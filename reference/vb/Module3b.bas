Attribute VB_Name = "Module3"


Public Sub Build3SeqTable()

Form1.SSPanel1 = "Building 3Seq P-value lookup table"
'Find max memory needed by looking at max permdist()
Dim MaxPermDists As Long, w As Long, X As Long, Y As Long, z As Long, YTable() As Single
ReDim WTable(1, 1, 1, 1)
MaxPermDists = 0
For X = 0 To Nextno - 1
    For Y = X + 1 To Nextno
        If PermDiffs(X, Y) > MaxPermDists Then
            MaxPermDists = PermDiffs(X, Y)
        End If
    Next Y
Next X


'Work out how much memory is available
Dim MaxTableSize As Double, SizeStore As Long
GlobalMemoryStatus MemSit
MaxTableSize = MemSit.dwTotalPhys


'MaxTableSize = MaxTableSize / 2 'only allow use of 1/2 of available memory

MaxTableSize = MaxTableSize / 16 '4 bytes for single,only allow use of 1/4 of phys mem
'XX = 1000 ^ (1 / 4)
MaxTableSize = MaxTableSize ^ (1 / 4) '4dimentions
MaxTableSize = MaxTableSize / 2 'only permit use of 1/2 of memory


If MaxPermDists > MaxTableSize Then
    MaxPermDists = CLng(MaxTableSize)
End If

ODir = CurDir
ChDir App.Path

ReDim XTable(0, 0, 0)

Open "3seqTable" For Binary As #1
Get #1, , SizeStore
If MaxPermDists <= SizeStore Then
    ReDim XTable(SizeStore + 1, SizeStore + 1, SizeStore + 1)
    Get #1, , XTable()
End If
Close #1
ChDir ODir

If UBound(XTable, 1) = 0 Then 'calculate the table if a big enough table is not available
    
StepBackin:


    On Error GoTo DropMaxPermDists
    ReDim XTable(MaxPermDists + 1, MaxPermDists + 1, MaxPermDists + 1)
    Form1.ProgressBar1 = 3
    ReDim YTable(MaxPermDists + 1, MaxPermDists + 1, MaxPermDists + 1, MaxPermDists + 1)
    
    On Error GoTo 0
    Form1.ProgressBar1 = 10
    Dummy = Fill4DSingArray(MaxPermDists + 1, MaxPermDists + 1, MaxPermDists + 1, MaxPermDists + 1, -1, YTable(0, 0, 0, 0))
    Form1.ProgressBar1 = 15
    Tot = (MaxPermDists + 1) ^ 3
    
    
    
    'AbortFlag = 1
    X = 0
    For nM = 0 To MaxPermDists
        For nN = 0 To MaxPermDists
            For nK = 0 To MaxPermDists
                dPValue = Seq3PVals(UBound(YTable, 1), UBound(YTable, 3), nM, nN, nK, YTable(0, 0, 0, 0))
                XTable(nM, nN, nK) = dPValue
                X = X + 1
            Next
        Next
       
        Form1.ProgressBar1 = 15 + X / Tot * 85
    Next
    ODir = CurDir
    ChDir App.Path
    Open "3seqTable" For Binary As #1
    Put #1, , MaxPermDists
    Put #1, , XTable()
    Close #1
    ChDir ODir
End If

XTableFlag = 1
Form1.ProgressBar1 = 0
Form1.SSPanel1 = ""
Exit Sub
DropMaxPermDists:
MaxPermDists = CLng(MaxPermDists / 2)
ReDim XTable(0, 0, 0), YTable(0, 0, 0, 0)
On Error GoTo 0
GoTo StepBackin
End Sub

Public Sub GetTSPVal(WF, nM, nN, nK, dPValue)
Dim PVM As Double, onM As Double, onN As Double
onN = nN
onM = nM
onk = nK
If nM >= UBound(XTable, 1) - 1 Or nN >= UBound(XTable, 1) - 1 Or nK >= UBound(XTable, 1) - 1 Then
    If nM >= nN And nM >= nK Then
        PVM = nM / (UBound(XTable, 1) - 2)
    ElseIf nN >= nM And nN >= nK Then
        PVM = nN / (UBound(XTable, 1) - 2)
    Else
        PVM = nK / (UBound(XTable, 1) - 2)
    End If
    
    nM = Int(nM / PVM)
    nN = Int(nN / PVM)
    nK = Int(nK / PVM)
    If nM > 0 Then
        PVM = onM / nM
        If onM / nM > PVM Then
            PVM = onM / nM
        End If
    ElseIf nN > 0 Then
        PVM = onN / nN
        If onN / nN > PVM Then
            PVM = onN / nN
        End If
    End If
    
    WF = 1
Else
    WF = 0
    PVM = 1
End If

dPValue = XTable(nM, nN, nK)
odPValue = dPValue
If PVM > 1 Then
    dPValue = dPValue ^ PVM ' / (10 ^ (PVM - 1))
End If
If dPValue = 0 And odPValue > 0 Then
    dPValue = 10 ^ -300
End If
nM = onM
nN = onN
nK = onk
End Sub
Public Sub WriteNames2(tseq1, tseq2, tseq3, red, green, blue)
OFontSize = Form1.Picture7.FontSize
L1 = Form1.Picture7.TextWidth(StraiName(tseq1))
L2 = Form1.Picture7.TextWidth(StraiName(tseq2))
L3 = Form1.Picture7.TextWidth(StraiName(tseq3))
DrawLen = Form1.Picture7.ScaleWidth - 30
TotLen = L1 + L2 + L3

If DrawLen > TotLen Then
    LOSpace = (DrawLen - TotLen) / 2
    XPos1 = 25
    XPos2 = XPos1 + L1 + LOSpace
    XPos3 = XPos2 + L2 + LOSpace
Else

    Do Until TotLen < DrawLen
        Form1.Picture7.FontSize = Form1.Picture7.FontSize - 1
        L1 = Form1.Picture7.TextWidth(StraiName(tseq1))
        L2 = Form1.Picture7.TextWidth(StraiName(tseq2))
        L3 = Form1.Picture7.TextWidth(StraiName(tseq3))
        TotLen = L1 + L2 + L3
    Loop

    LOSpace = (DrawLen - TotLen) / 2
    XPos1 = 25
    XPos2 = XPos1 + L1 + LOSpace
    XPos3 = XPos2 + L2 + LOSpace
End If
YPos = Form1.Picture7.Height * (0.92)
Form1.Picture7.ForeColor = ThreeQuaterColour
Form1.Picture7.CurrentX = XPos1 - 1
Form1.Picture7.CurrentY = YPos - 1
Form1.Picture7.Print StraiName(tseq1)
Form1.Picture7.CurrentX = XPos2 - 1
Form1.Picture7.CurrentY = YPos - 1
Form1.Picture7.Print StraiName(tseq2)
Form1.Picture7.CurrentX = XPos3 - 1
Form1.Picture7.CurrentY = YPos - 1
Form1.Picture7.Print StraiName(tseq3)
Form1.Picture7.ForeColor = QuaterColour
Form1.Picture7.CurrentX = XPos1 + 1
Form1.Picture7.CurrentY = YPos + 1
Form1.Picture7.Print StraiName(tseq1)
Form1.Picture7.CurrentX = XPos2 + 1
Form1.Picture7.CurrentY = YPos + 1
Form1.Picture7.Print StraiName(tseq2)
Form1.Picture7.CurrentX = XPos3 + 1
Form1.Picture7.CurrentY = YPos + 1
Form1.Picture7.Print StraiName(tseq3)
Form1.Picture7.CurrentX = XPos1
Form1.Picture7.CurrentY = YPos
Form1.Picture7.ForeColor = red 'RGB(255, 0, 0)'this is actually green
Form1.Picture7.Print StraiName(tseq1)
Form1.Picture7.CurrentX = XPos2
Form1.Picture7.CurrentY = YPos
Form1.Picture7.ForeColor = blue 'RGB(0, 255, 0)'this is actually blue
Form1.Picture7.Print StraiName(tseq2)
Form1.Picture7.CurrentX = XPos3
Form1.Picture7.CurrentY = YPos
Form1.Picture7.ForeColor = green 'RGB(0, 0, 255)'this is actually red
Form1.Picture7.Print StraiName(tseq3)
'Form1.Picture7.Refresh
Form1.Picture7.FontSize = OFontSize
End Sub
Public Sub TSXOverC()
Dim BE As Long, EN As Long, BE2 As Long, EN2 As Long, nK As Long, nM As Long, nN As Long, nL As Long, X As Long, CurrentHeight As Long, XoverSeqNumTS() As Long, XOverSeqNumPerms() As Long, Y As Long, MaxSeen As Long, MaxDescentSeen As Long
Dim ii As Long, dPValue As Double
tseq1 = Seq1
tseq2 = Seq2
tseq3 = Seq3

If XTableFlag = 0 Then Call Build3SeqTable

Seq1 = tseq3 ' this silly conversion needed so names and colours match up with those of chimaera
Seq2 = tseq2
Seq3 = tseq1


'Call TSXOver(0)
ReDim XoverSeqNumTS(Len(StrainSeq(0)))


Dim tXOverSeqNum() As Long

Rnd (GetTickCount)
Dim Tmp As Long, THold(2) As Long, SeqRnd() As Integer, dPValue2 As Double, xPosdiffx() As Long, xDiffposx() As Long, sn As Long, sm As Long, sk As Long, sbe As Long, sen As Long


        
Dim BES1 As Long, BES2 As Long, BES3 As Long, ENS1 As Long, ENS2 As Long, ENS3 As Long
Dim BE2S1 As Long, BE2S2 As Long, BE2S3 As Long, EN2S1 As Long, EN2S2 As Long, EN2S3 As Long
Dim nMS1 As Long, nNS1 As Long, nKS1 As Long, nLS1 As Long, nMS2 As Long, nNS2 As Long, nKS2 As Long, nLS2 As Long, nMS3 As Long, nNS3 As Long, nKS3 As Long, nLS3 As Long
Dim XPosDiffS1() As Long, XDiffposS1() As Long, XPosDiffS2() As Long, XDiffposS2() As Long, XPosDiffS3() As Long, XDiffposS3() As Long
Dim XoverSeqNumTS1() As Long, XoverSeqNumTS2() As Long, XoverSeqNumTS3() As Long
Dim LenXoverSeqS1 As Long, LenXoverSeqS2 As Long, LenXoverSeqS3 As Long

ReDim XPosDiffS1(Len(StrainSeq(0))), XDiffposS1(Len(StrainSeq(0))), XPosDiffS2(Len(StrainSeq(0))), XDiffposS2(Len(StrainSeq(0))), XPosDiffS3(Len(StrainSeq(0))), XDiffposS3(Len(StrainSeq(0)))
ReDim XoverSeqNumTS1(Len(StrainSeq(0))), XoverSeqNumTS2(Len(StrainSeq(0))), XoverSeqNumTS3(Len(StrainSeq(0)))

LenXoverSeqS1 = FindSubSeqTS2(Len(StrainSeq(0)), Seq2, Seq3, Seq1, BES1, ENS1, BE2S1, EN2S1, nMS1, nNS1, nKS1, nLS1, XPosDiffS1(0), XDiffposS1(0), SeqNum(0, 0), XoverSeqNumTS1(0), MissingData(0, 0))
LenXoverSeqS2 = FindSubSeqTS2(Len(StrainSeq(0)), Seq1, Seq3, Seq2, BES2, ENS2, BE2S2, EN2S2, nMS2, nNS2, nKS2, nLS2, XPosDiffS2(0), XDiffposS2(0), SeqNum(0, 0), XoverSeqNumTS2(0), MissingData(0, 0))
LenXoverSeqS3 = FindSubSeqTS2(Len(StrainSeq(0)), Seq1, Seq2, Seq3, BES3, ENS3, BE2S3, EN2S3, nMS3, nNS3, nKS3, nLS3, XPosDiffS3(0), XDiffposS3(0), SeqNum(0, 0), XoverSeqNumTS3(0), MissingData(0, 0))
'2,7,3:2642,1208,1208,2642:26,20,12,19
        
        '1312,2309
        '2352,1306
        
        '1306,2309
        '2309,1306
'45 26,20,12,19
If X = X Then
    Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeqS1, 1, CircularFlag, nKS1, BES1, ENS1, XDiffposS1(0), XPosDiffS1(0), XoverSeqNumTS1(0))
    Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeqS1, -1, CircularFlag, nLS1, BE2S1, EN2S1, XDiffposS1(0), XPosDiffS1(0), XoverSeqNumTS1(0))

    Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeqS2, 1, CircularFlag, nKS2, BES2, ENS2, XDiffposS2(0), XPosDiffS2(0), XoverSeqNumTS2(0))
    Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeqS2, -1, CircularFlag, nLS2, BE2S2, EN2S2, XDiffposS2(0), XPosDiffS2(0), XoverSeqNumTS2(0))

    Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeqS3, 1, CircularFlag, nKS3, BES3, ENS3, XDiffposS3(0), XPosDiffS3(0), XoverSeqNumTS3(0))
    Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeqS3, -1, CircularFlag, nLS3, BE2S3, EN2S3, XDiffposS3(0), XPosDiffS3(0), XoverSeqNumTS3(0))
    '12,2642,1208------19,2653,1208
    '19,1208,2642------19,1214,2642
    X = X
Else
    Call CheckWrap(LenXoverSeqS1, nKS1, BES1, ENS1, 1, XDiffposS1(), XPosDiffS1(), XoverSeqNumTS1())
    Call CheckWrap(LenXoverSeqS1, nLS1, BE2S1, EN2S1, -1, XDiffposS1(), XPosDiffS1(), XoverSeqNumTS1())
            
    Call CheckWrap(LenXoverSeqS2, nKS2, BES2, ENS2, 1, XDiffposS2(), XPosDiffS2(), XoverSeqNumTS2())
    Call CheckWrap(LenXoverSeqS2, nLS2, BE2S2, EN2S2, -1, XDiffposS2(), XPosDiffS2(), XoverSeqNumTS2())
    
    Call CheckWrap(LenXoverSeqS3, nKS3, BES3, ENS3, 1, XDiffposS3(), XPosDiffS3(), XoverSeqNumTS3())
    Call CheckWrap(LenXoverSeqS3, nLS3, BE2S3, EN2S3, -1, XDiffposS3(), XPosDiffS3(), XoverSeqNumTS3())
End If
'13,2653,1208
'19,1214,2642
       
Dim dPValueS1 As Double, dPValueS2 As Double, dPValueS3 As Double
Dim dPValue2S1 As Double, dPValue2S2 As Double, dPValue2S3 As Double
Dim WFS1, WF2S1, WFS2, WF2S2, WFS3, WF2S3

GetTSPVal WFS1, nMS1, nNS1, nKS1, dPValueS1
GetTSPVal WF2S1, nNS1, nMS1, nLS1, dPValue2S1

GetTSPVal WFS2, nMS2, nNS2, nKS2, dPValueS2
GetTSPVal WF2S2, nNS2, nMS2, nLS2, dPValue2S2

GetTSPVal WFS3, nMS3, nNS3, nKS3, dPValueS3
GetTSPVal WF2S3, nNS3, nMS3, nLS3, dPValue2S3
        
Dim MultNegS1, MultNegS2, MultNegS3
        
MultNegS1 = 0: MultNegS2 = 0: MultNegS3 = 0

If dPValueS1 = 0 Then dPValueS1 = 10
If dPValueS2 = 0 Then dPValueS2 = 10
If dPValueS3 = 0 Then dPValueS3 = 10
If dPValue2S1 = 0 Then dPValue2S1 = 10
If dPValue2S2 = 0 Then dPValue2S2 = 10
If dPValue2S3 = 0 Then dPValue2S3 = 10


If dPValue2S1 < dPValueS1 Then
    Call SwapRound(WFS1, WF2S1, BES1, ENS1, BE2S1, EN2S1, nMS1, nNS1, nKS1, nLS1, dPValueS1, dPValue2S1, MultNegS1)
End If

If dPValue2S2 < dPValueS2 Then
    Call SwapRound(WFS2, WF2S2, BES2, ENS2, BE2S2, EN2S2, nMS2, nNS2, nKS2, nLS2, dPValueS2, dPValue2S2, MultNegS2)
End If

If dPValue2S3 < dPValueS3 Then
    Call SwapRound(WF3, WF2S3, BES3, ENS3, BE2S3, EN2S3, nMS3, nNS3, nKS3, nLS3, dPValueS3, dPValue2S3, MultNegS3)
End If


If MultNegS1 = 1 Then
    For X = 0 To LenXoverSeqS1
        XoverSeqNumTS1(X) = XoverSeqNumTS1(X) * -1
    Next X
End If

If MultNegS2 = 1 Then
    For X = 0 To LenXoverSeqS2
        XoverSeqNumTS2(X) = XoverSeqNumTS2(X) * -1
    Next X
End If

If MultNegS3 = 1 Then
    For X = 0 To LenXoverSeqS3
        XoverSeqNumTS3(X) = XoverSeqNumTS3(X) * -1
    Next X
End If

    '1.646-4,699,18
    
Call CheckSplit3Seq(BES1, ENS1, nMS1, nNS1, dPValueS1, LenXoverSeqS1, XoverSeqNumTS1(), XPosDiffS1(), XDiffposS1())

If dPValueS1 > dPValue2S1 Or FindallFlag = 1 Then
    
    If MultNegS1 = 0 Then
        For X = 0 To LenXoverSeqS1
            XoverSeqNumTS1(X) = XoverSeqNumTS1(X) * -1
        Next X
        'Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
    End If
    'Call SwapRound(WFS1, WF2S1, BES1, ENS1, BE2S1, EN2S1, nMS1, nNS1, nKS1, nLS1, dPValueS1, dPValue2S1, MultNegS1)
    
    
    
    Call CheckSplit3Seq(BE2S1, EN2S1, nNS1, nMS1, dPValue2S1, LenXoverSeqS1, XoverSeqNumTS1(), XPosDiffS1(), XDiffposS1())
    If FindallFlag = 0 And dPValueS1 > dPValue2S1 Then
        If MultNegS1 = 1 Then
            MultNegS1 = 0
            
        Else
            MultNegS1 = 1
            
            
        End If
        
        Call SwapRound(WFS1, WF2S1, BES1, ENS1, BE2S1, EN2S1, nMS1, nNS1, nKS1, nLS1, dPValueS1, dPValue2S1, MultNegS1)
    End If

End If
        
        
Call CheckSplit3Seq(BES2, ENS2, nMS2, nNS2, dPValueS2, LenXoverSeqS2, XoverSeqNumTS2(), XPosDiffS2(), XDiffposS2())

If dPValueS2 > dPValue2S2 Or FindallFlag = 1 Then
    If MultNegS2 = 0 Then
        For X = 0 To LenXoverSeqS2
            XoverSeqNumTS2(X) = XoverSeqNumTS2(X) * -1
        Next X
        'Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
    End If
    
    Call CheckSplit3Seq(BE2S2, EN2S2, nNS2, nMS2, dPValue2S2, LenXoverSeqS2, XoverSeqNumTS2(), XPosDiffS2(), XDiffposS2())
    If FindallFlag = 0 And dPValueS2 > dPValue2S2 Then
        If MultNegS2 = 1 Then
            MultNegS2 = 0
            
        Else
            MultNegS2 = 1
        End If
        Call SwapRound(WFS2, WF2S2, BES2, ENS2, BE2S2, EN2S2, nMS2, nNS2, nKS2, nLS2, dPValueS2, dPValue2S2, MultNegS2)
        
    End If

End If

Call CheckSplit3Seq(BES3, ENS3, nMS3, nNS3, dPValueS3, LenXoverSeqS3, XoverSeqNumTS3(), XPosDiffS3(), XDiffposS3())

If dPValueS3 > dPValue2S3 Or FindallFlag = 1 Then
    If MultNegS3 = 0 Then
        For X = 0 To LenXoverSeqS3
            XoverSeqNumTS3(X) = XoverSeqNumTS3(X) * -1
        Next X
        'Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
    End If
    
    Call CheckSplit3Seq(BE2S3, EN2S3, nNS3, nMS3, dPValue2S3, LenXoverSeqS3, XoverSeqNumTS3(), XPosDiffS3(), XDiffposS3())
    If FindallFlag = 0 And dPValueS3 > dPValue2S3 Then
        If MultNegS3 = 1 Then
            MultNegS3 = 0
            
        Else
            MultNegS3 = 1
            
            
        End If
        Call SwapRound(WFS3, WF2S3, BES3, ENS3, BE2S3, EN2S3, nMS3, nNS3, nKS3, nLS3, dPValueS3, dPValue2S3, MultNegS3)
        
    End If

End If
        
        
        
If BES1 = 0 Then BES1 = XDiffposS1(0)
If BES2 = 0 Then BES2 = XDiffposS2(0)
If BES3 = 0 Then BES3 = XDiffposS3(0)

If LenXoverSeqS1 < 1 Or LenXoverSeqS2 < 1 Or LenXoverSeqS2 < 1 Then Exit Sub


Dim tMissingData() As Long
ReDim tMissingData(Len(StrainSeq(0)), 2)
For X = 0 To Len(StrainSeq(0))
    tMissingData(X, 0) = MissingData(X, Seq1)
    tMissingData(X, 1) = MissingData(X, Seq2)
    tMissingData(X, 2) = MissingData(X, Seq3)
Next X

Dim PermRangeS1() As Long, PermRangeS2() As Long, PermRangeS3() As Long
ReDim PermRangeS1(1, LenXoverSeqS1), PermRangeS2(1, LenXoverSeqS2), PermRangeS3(1, LenXoverSeqS3)

For A = 0 To 2
    ReDim XOverSeqNumPerms(Len(StrainSeq(0)), 100)
    ReDim tXOverSeqNum(Len(StrainSeq(0)))
    ReDim SeqRnd(Len(StrainSeq(0)), 2)
    ReDim PermRange(1, Len(StrainSeq(0)))
    
    If A = 0 Then
        If MultNegS1 = 0 Then
            Seq1 = tseq2
            Seq2 = tseq1
        Else
            Seq2 = tseq2
            Seq1 = tseq1
        End If
        Seq3 = tseq3
    ElseIf A = 1 Then
        
        If MultNegS2 = 0 Then
            Seq1 = tseq3
            Seq2 = tseq1
        Else
            Seq2 = tseq3
            Seq1 = tseq1
        End If
        Seq3 = tseq2
    Else
        If MultNegS3 = 0 Then
            Seq1 = tseq3
            Seq2 = tseq2
        Else
            Seq2 = tseq3
            Seq1 = tseq2
        End If
        
        Seq3 = tseq1
    End If
    For X = 1 To 100
        ReDim xPosdiffx(Len(StrainSeq(0))), xDiffposx(Len(StrainSeq(0)))
        For z = 1 To Len(StrainSeq(0))
            
            SeqRnd(z, 0) = SeqNum(z, Seq1)
            SeqRnd(z, 1) = SeqNum(z, Seq2)
            SeqRnd(z, 2) = SeqNum(z, Seq3)
        Next z
       
        For z = 1 To Len(StrainSeq(0))
            NewPos = CLng(Rnd * (Len(StrainSeq(0)) - 1)) + 1
            If tMissingData(z, 0) = 0 And tMissingData(z, 1) = 0 And tMissingData(z, 2) = 0 And tMissingData(NewPos, 0) = 0 And tMissingData(NewPos, 1) = 0 And tMissingData(NewPos, 2) = 0 Then
                THold(0) = SeqRnd(z, 0)
                THold(1) = SeqRnd(z, 1)
                THold(2) = SeqRnd(z, 2)
                SeqRnd(z, 0) = SeqRnd(NewPos, 0)
                SeqRnd(z, 1) = SeqRnd(NewPos, 1)
                SeqRnd(z, 2) = SeqRnd(NewPos, 2)
                SeqRnd(NewPos, 0) = THold(0)
                SeqRnd(NewPos, 1) = THold(1)
                SeqRnd(NewPos, 2) = THold(2)
            End If
        Next z
        ReDim tXOverSeqNum(Len(StrainSeq(0)))
        sbe = 0
        sen = 0
        sm = 0
        sn = 0
        sk = 0
        
        LenXoverSeqx = FindSubSeqTS2(Len(StrainSeq(0)), 0, 1, 2, sbe, sen, BE2, EN2, sm, sn, sk, nL, xPosdiffx(0), xDiffposx(0), SeqRnd(0, 0), tXOverSeqNum(0), tMissingData(0, 0))
        For z = 0 To Len(StrainSeq(0))
            XOverSeqNumPerms(z, X) = tXOverSeqNum(z)
        Next z
    Next X
     
    For X = 0 To LenXoverSeqx
        PermRange(1, X) = LenXoverSeqx
        PermRange(0, X) = -LenXoverSeqx
        For Y = 1 To 100
            If XOverSeqNumPerms(X, Y) > PermRange(0, X) Then PermRange(0, X) = XOverSeqNumPerms(X, Y)
            If XOverSeqNumPerms(X, Y) < PermRange(1, X) Then PermRange(1, X) = XOverSeqNumPerms(X, Y)
        Next Y
    Next X
    If A = 0 Then
        For X = 0 To LenXoverSeqS1
            PermRangeS1(1, X) = PermRange(1, X)
            PermRangeS1(0, X) = PermRange(0, X)
        Next X
    ElseIf A = 1 Then
        For X = 0 To LenXoverSeqS2
            PermRangeS2(1, X) = PermRange(1, X)
            PermRangeS2(0, X) = PermRange(0, X)
        Next X
    Else
        For X = 0 To LenXoverSeqS3
            PermRangeS3(1, X) = PermRange(1, X)
            PermRangeS3(0, X) = PermRange(0, X)
        Next X
    End If
Next A
    
 If dPValueS1 > 0 And (dPValueS1 <= dPValueS2 Or dPValueS2 = 0) And (dPValueS1 <= dPValueS3 Or dPValueS3 = 0) Then
    dPValue = dPValueS1
    EN = ENS1
    BE = BES1
 ElseIf dPValueS2 > 0 And (dPValueS2 <= dPValueS1 Or dPValueS1 = 0) And (dPValueS2 <= dPValueS3 Or dPValueS3 = 0) Then
    dPValue = dPValueS2
    EN = ENS2
    BE = BES2
 Else
    dPValue = dPValueS3
    EN = ENS3
    BE = BES3
 End If
' performs a Dunn-Sidak correction for pval with m trials
 xpvalue = dPValue * MCCorrection
 If dPValue >= 1 Then
     dPValue = 1
 Else
     dPValue = 1 - (1 - dPValue) ^ MCCorrection
 End If
 'characterise the event
 If dPValue = 0 Then dPValue = xpvalue
 Dim Max As Double, Min As Double
 Max = -LenXoverSeq
 Min = LenXoverSeq
 For X = 0 To LenXoverSeqS1
     If Max < XoverSeqNumTS1(X) Then Max = XoverSeqNumTS1(X)
     If Min > XoverSeqNumTS1(X) Then Min = XoverSeqNumTS1(X)
 Next X
 
For X = 0 To LenXoverSeqS2
    If Max < XoverSeqNumTS2(X) Then Max = XoverSeqNumTS2(X)
    If Min > XoverSeqNumTS2(X) Then Min = XoverSeqNumTS2(X)
Next X

For X = 0 To LenXoverSeqS3
    If Max < XoverSeqNumTS3(X) Then Max = XoverSeqNumTS3(X)
    If Min > XoverSeqNumTS3(X) Then Min = XoverSeqNumTS3(X)
Next X


For X = 0 To LenXoverSeqS1
    If Max < PermRangeS1(0, X) Then Max = PermRangeS1(0, X)
    If Min > PermRangeS1(1, X) Then Min = PermRangeS1(1, X)
Next X

For X = 0 To LenXoverSeqS2
    If Max < PermRangeS2(0, X) Then Max = PermRangeS2(0, X)
    If Min > PermRangeS2(1, X) Then Min = PermRangeS2(1, X)
Next X

For X = 0 To LenXoverSeqS3
    If Max < PermRangeS3(0, X) Then Max = PermRangeS3(0, X)
    If Min > PermRangeS3(1, X) Then Min = PermRangeS3(1, X)
Next X
 
 Dim YScaleFactor As Double, PntAPI As POINTAPI, Pict As Long
Form1.Picture7.Picture = LoadPicture()
YScaleFactor = 0.85
PicHeight = Form1.Picture7.Height * YScaleFactor
XFactor = ((Form1.Picture7.Width - 40) / Len(StrainSeq(0)))


 Dim red As Long, green As Long, blue As Long, hRed As Long, hGreen As Long, hBlue As Long


Call DoAxes(0, Len(StrainSeq(0)), -1, Max, Min, 0, "Height")




red = RGB(CInt(BkR - BkR / 2), CInt(BkG - BkG / 6), CInt(BkB - BkB / 2))
blue = RGB(CInt(BkR - BkR / 6), CInt(BkG - BkG / 2), CInt(BkB - BkB / 2))
green = RGB(CInt(BkR - BkR / 2), CInt(BkG - BkG / 2), CInt(BkB - BkB / 6))


If TManFlag <> 22 Then
    Call Highlight
Else
    BE = 0
    EN = Len(StrainSeq(0))
    GBlockNum = -1
End If

'Colour the bkround white



hRed = RGB(254, 254, 254)
hBlue = RGB(253, 253, 253)
hGreen = RGB(252, 252, 252)
Form1.Picture7.DrawMode = 13
Form1.Picture7.ForeColor = hGreen + 1
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + XDiffposS1(0) * XFactor, (PicHeight - (15 + ((PermRangeS1(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
For X = 0 To LenXoverSeqS1
    LineTo Pict, 30 + XDiffposS1(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS1(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X
For X = LenXoverSeqS1 To 0 Step -1
    LineTo Pict, 30 + XDiffposS1(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS1(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X

LineTo Pict, 30 + XDiffposS1(0) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS1(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))

Form1.Picture7.FillStyle = 0
Form1.Picture7.FillColor = hGreen
Pict = Form1.Picture7.hDC
'dy = PicHeight - (15 + ((((PermRangeS1(1, 0) + (PermRangeS1(0, 0) - PermRangeS1(1, 0)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
'FloodFill Pict, 30 + XDiffposS1(0) * XFactor + XFactor + 2, dy, hGreen + 1
X = 0
GoOn = 0

Dim MaxR As Long, WinR As Long
WinR = 10000000
MaxR = 0
Do While X <= LenXoverSeqS1
    If PermRangeS1(0, X) - PermRangeS1(1, X) > MaxR Then
        MaxR = PermRangeS1(0, X) - PermRangeS1(1, X)
        WinR = X
    End If
    X = X + 1
Loop

X = WinR


'Do While X <= LenXoverSeqS1
    If PermRangeS1(1, X) <> PermRangeS1(0, X) Then
        If X < LenXoverSeqS1 Then
            dy = PicHeight - (15 + ((((PermRangeS1(1, X) + (PermRangeS1(0, X) - PermRangeS1(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
            dx = 30 + XDiffposS1(X) * XFactor + XFactor
        Else
            dy = PicHeight - (15 + ((((PermRangeS1(1, X) + (PermRangeS1(0, X) - PermRangeS1(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
            dx = 30 + XDiffposS1(X) * XFactor + XFactor
        End If
        GoOn = 1
    End If
'    X = X + 1
'Loop
'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
Pict = Form1.Picture7.hDC

If GoOn = 1 Then
    FloodFill Pict, dx, dy, hGreen + 1
End If

Form1.Picture7.FillStyle = 1


Form1.Picture7.ForeColor = hBlue + 1
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + XDiffposS2(0) * XFactor, (PicHeight - (15 + ((PermRangeS2(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
For X = 0 To LenXoverSeqS2
    LineTo Pict, 30 + XDiffposS2(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS2(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X
For X = LenXoverSeqS2 To 0 Step -1
    LineTo Pict, 30 + XDiffposS2(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS2(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X

LineTo Pict, 30 + XDiffposS2(0) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS2(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))

Form1.Picture7.FillStyle = 0
Form1.Picture7.FillColor = hBlue
Pict = Form1.Picture7.hDC
'dy = PicHeight - (15 + ((((PermRangeS2(1, 0) + (PermRangeS2(0, 0) - PermRangeS2(1, 0)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
'FloodFill Pict, 30 + XDiffposS2(0) * XFactor + XFactor + 2, dy, hBlue + 1
Form1.Picture7.FillStyle = 1
X = 0
GoOn = 0
WinR = 10000000
MaxR = 0
Do While X <= LenXoverSeqS2
    If PermRangeS2(0, X) - PermRangeS2(1, X) > MaxR Then
        MaxR = PermRangeS2(0, X) - PermRangeS2(1, X)
        WinR = X
    End If
    X = X + 1
Loop

X = WinR

'Do While X <= LenXoverSeqS2
    If PermRangeS2(1, X) <> PermRangeS2(0, X) Then
        If X < LenXoverSeqS2 Then
            dy = PicHeight - (15 + ((((PermRangeS2(1, X) + (PermRangeS2(0, X) - PermRangeS2(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
            dx = 30 + XDiffposS2(X) * XFactor + XFactor
        Else
            dy = PicHeight - (15 + ((((PermRangeS2(1, X) + (PermRangeS2(0, X) - PermRangeS2(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
            dx = 30 + XDiffposS2(X) * XFactor + XFactor
        End If
        GoOn = 1
    End If
'    X = X + 1
'Loop
'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
'Form1.Picture7.AutoRedraw = True
Form1.Picture7.FillStyle = 0
Form1.Picture7.FillColor = hBlue
Pict = Form1.Picture7.hDC
If GoOn = 1 Then
    FloodFill Pict, CLng(dx), CLng(dy), hBlue + 1
End If
'Form1.Picture7.Circle (dx, dy), 10, 0
'Form1.Picture7.ForeColor = hBlue + 1
'Form1.Picture7.DrawWidth = 1

'MoveToEx Pict, dx, dy, PntAPI
'LineTo Pict, dx, dy

Form1.Picture7.ForeColor = hRed + 1
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + XDiffposS3(0) * XFactor, (PicHeight - (15 + ((PermRangeS3(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
For X = 0 To LenXoverSeqS3
    LineTo Pict, 30 + XDiffposS3(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS3(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X
For X = LenXoverSeqS3 To 0 Step -1
    LineTo Pict, 30 + XDiffposS3(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS3(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X

LineTo Pict, 30 + XDiffposS3(0) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS3(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))
Form1.Picture7.FillStyle = 0
Form1.Picture7.FillColor = hRed
Pict = Form1.Picture7.hDC
'dy = PicHeight - (15 + ((((PermRangeS3(1, 0) + (PermRangeS3(0, 0) - PermRangeS3(1, 0)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
'FloodFill Pict, 30 + XDiffposS3(0) * XFactor + XFactor + 2, dy, hRed + 1

X = 0
WinR = 10000000
MaxR = 0
Do While X <= LenXoverSeqS3
    If PermRangeS3(0, X) - PermRangeS3(1, X) > MaxR Then
        MaxR = PermRangeS3(0, X) - PermRangeS3(1, X)
        WinR = X
    End If
    X = X + 1
Loop

X = WinR
GoOn = 0
'Do While X <= LenXoverSeqS3
    If PermRangeS3(1, X) <> PermRangeS3(0, X) Then
        If X < LenXoverSeqS3 Then
            dy = PicHeight - (15 + ((((PermRangeS3(1, X) + (PermRangeS3(0, X) - PermRangeS3(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
            dx = 30 + XDiffposS3(X) * XFactor + XFactor
        Else
            dy = PicHeight - (15 + ((((PermRangeS3(1, X) + (PermRangeS3(0, X) - PermRangeS3(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
            dx = 30 + XDiffposS3(X) * XFactor + XFactor
        End If
        GoOn = 1
    End If
'    X = X + 1
'Loop
'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
Pict = Form1.Picture7.hDC

If GoOn = 1 Then
    FloodFill Pict, dx, dy, hRed + 1
End If

XX = Form1.Picture7.FillStyle


Form1.Picture7.FillStyle = 1
Form1.Picture7.DrawMode = 13



'Now do the plots
'hRed = RGB(CInt(BkR - BkR / 4), CInt(BkG - BkG / 24), CInt(BkB - BkB / 4)) '
'    hBlue = RGB(CInt(BkR - BkR / 24), CInt(BkB - BkB / 4), CInt(BkG - BkG / 4))
'    hGreen = RGB(CInt(BkR - BkR / 4), CInt(BkG - BkG / 4), CInt(BkB - BkB / 24))

hRed = RGB(CInt(BkR - BkR / 4), CInt(BkG - BkG / 24), CInt(BkB - BkB / 3))
hBlue = RGB(CInt(BkR - BkR / 24), CInt(BkB - BkB / 3), CInt(BkG - BkG / 4))
hGreen = RGB(CInt(BkR - BkR / 3), CInt(BkG - BkG / 4), CInt(BkB - BkB / 24))

For z = 0 To 0
    If z = 0 Then
        curm = 9
    ElseIf z = 1 Then
        
        curm = 14
    ElseIf z = 2 Then
        curm = 10
    End If
     '6,14,10
     '3,6,9,5,5
    
    Form1.Picture7.FillStyle = 1
    Form1.Picture7.DrawMode = 13
    Form1.Picture7.ForeColor = hRed + 1
    Pict = Form1.Picture7.hDC
    MoveToEx Pict, 30 + XDiffposS3(0) * XFactor, (PicHeight - (15 + ((PermRangeS3(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
    For X = 0 To LenXoverSeqS3
        LineTo Pict, 30 + XDiffposS3(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS3(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    For X = LenXoverSeqS3 To 0 Step -1
        LineTo Pict, 30 + XDiffposS3(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS3(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    
    LineTo Pict, 30 + XDiffposS3(0) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS3(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))
    Form1.Picture7.DrawMode = curm
    Form1.Picture7.FillStyle = 0
    Form1.Picture7.FillColor = hRed
    Pict = Form1.Picture7.hDC
    
    'dy = PicHeight - (15 + ((((PermRangeS3(1, 0) + (PermRangeS3(0, 0) - PermRangeS3(1, 0)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
    
    X = 0
    WinR = 10000000
    MaxR = 0
    Do While X <= LenXoverSeqS3
        If PermRangeS3(0, X) - PermRangeS3(1, X) > MaxR Then
            MaxR = PermRangeS3(0, X) - PermRangeS3(1, X)
            WinR = X
        End If
        X = X + 1
    Loop
    
    X = WinR
    GoOn = 0
    'Do While X <= LenXoverSeqS3
        If PermRangeS3(1, X) <> PermRangeS3(0, X) Then
            If X < LenXoverSeqS3 Then
                dy = PicHeight - (15 + ((((PermRangeS3(1, X) + (PermRangeS3(0, X) - PermRangeS3(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
                dx = 30 + XDiffposS3(X) * XFactor + XFactor
            Else
                dy = PicHeight - (15 + ((((PermRangeS3(1, X) + (PermRangeS3(0, X) - PermRangeS3(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
                dx = 30 + XDiffposS3(X) * XFactor + XFactor
            End If
            GoOn = 1
        End If
    '    X = X + 1
    'Loop
    'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
    Pict = Form1.Picture7.hDC

    If GoOn = 1 Then
        FloodFill Pict, dx, dy, hRed + 1
    End If
    'FloodFill Pict, 30 + XDiffposS3(0) * XFactor + XFactor + 2, dy, hRed + 1
    Form1.Picture7.FillStyle = 1
    Form1.Picture7.DrawMode = 13
    
    Form1.Picture7.DrawMode = 13
    Form1.Picture7.ForeColor = hGreen + 1
    Pict = Form1.Picture7.hDC
    MoveToEx Pict, 30 + XDiffposS1(0) * XFactor, (PicHeight - (15 + ((PermRangeS1(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
    For X = 0 To LenXoverSeqS1
        LineTo Pict, 30 + XDiffposS1(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS1(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    For X = LenXoverSeqS1 To 0 Step -1
        LineTo Pict, 30 + XDiffposS1(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS1(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    
    LineTo Pict, 30 + XDiffposS1(0) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS1(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))
    Form1.Picture7.DrawMode = curm
    Form1.Picture7.FillStyle = 0
    Form1.Picture7.FillColor = hGreen
    Pict = Form1.Picture7.hDC
    'dy = PicHeight - (15 + ((((PermRangeS1(1, 0) + (PermRangeS1(0, 0) - PermRangeS1(1, 0)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
    'FloodFill Pict, 30 + XDiffposS1(0) * XFactor + XFactor + 2, dy, hGreen + 1
    
    X = 0
    GoOn = 0
    WinR = 10000000
    MaxR = 0
    Do While X <= LenXoverSeqS1
        If PermRangeS1(0, X) - PermRangeS1(1, X) > MaxR Then
            MaxR = PermRangeS1(0, X) - PermRangeS1(1, X)
            WinR = X
        End If
        X = X + 1
    Loop
    
    X = WinR
    'Do While X <= LenXoverSeqS1
        If PermRangeS1(1, X) <> PermRangeS1(0, X) Then
            If X < LenXoverSeqS1 Then
                dy = PicHeight - (15 + ((((PermRangeS1(1, X) + (PermRangeS1(0, X) - PermRangeS1(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
                dx = 30 + XDiffposS1(X) * XFactor + XFactor
            Else
                dy = PicHeight - (15 + ((((PermRangeS1(1, X) + (PermRangeS1(0, X) - PermRangeS1(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
                dx = 30 + XDiffposS1(X) * XFactor + XFactor
            End If
            GoOn = 1
        End If
    '    X = X + 1
    'Loop
    'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
    
    'Form1.Picture7.FillStyle = 0
    'Form1.Picture7.FillColor = hGreen
    Pict = Form1.Picture7.hDC
    If GoOn = 1 Then
        FloodFill Pict, dx, dy, hGreen + 1
    End If
    'Form1.Picture7.DrawWidth = 4
    'Form1.Picture7.DrawMode = 13
    'Form1.Picture7.ForeColor = 0
    'Pict = Form1.Picture7.hdc
    'MoveToEx Pict, dx, dy, PntAPI
    'LineTo Pict, dx, dy
    Form1.Picture7.FillStyle = 1
    
    Form1.Picture7.DrawMode = 13
    Form1.Picture7.ForeColor = hBlue + 1 'actually red
    Pict = Form1.Picture7.hDC
    MoveToEx Pict, 30 + XDiffposS2(0) * XFactor, (PicHeight - (15 + ((PermRangeS2(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
    For X = 0 To LenXoverSeqS2
        LineTo Pict, 30 + XDiffposS2(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS2(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    For X = LenXoverSeqS2 To 0 Step -1
        LineTo Pict, 30 + XDiffposS2(X) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS2(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    
    LineTo Pict, 30 + XDiffposS2(0) * XFactor + XFactor, PicHeight - (15 + (((PermRangeS2(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))
    Form1.Picture7.DrawMode = curm
    Form1.Picture7.FillStyle = 0
    Form1.Picture7.FillColor = hBlue
    Pict = Form1.Picture7.hDC
    'dy = PicHeight - (15 + ((((PermRangeS2(1, 0) + (PermRangeS2(0, 0) - PermRangeS2(1, 0)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
    'FloodFill Pict, 30 + XDiffposS2(0) * XFactor + XFactor + 2, dy, hBlue + 1
    X = 0
    GoOn = 0
    WinR = 10000000
    MaxR = 0
    Do While X <= LenXoverSeqS2
        If PermRangeS2(0, X) - PermRangeS2(1, X) > MaxR Then
            MaxR = PermRangeS2(0, X) - PermRangeS2(1, X)
            WinR = X
        End If
        X = X + 1
    Loop
    
    X = WinR
    
    'Do While X <= LenXoverSeqS2
        If PermRangeS2(1, X) <> PermRangeS2(0, X) Then
            If X < LenXoverSeqS2 Then
                dy = PicHeight - (15 + ((((PermRangeS2(1, X) + (PermRangeS2(0, X) - PermRangeS2(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
                dx = 30 + XDiffposS2(X) * XFactor + XFactor
            Else
                dy = PicHeight - (15 + ((((PermRangeS2(1, X) + (PermRangeS2(0, X) - PermRangeS2(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
                dx = 30 + XDiffposS2(X) * XFactor + XFactor
            End If
            GoOn = 1
        End If
    '    X = X + 1
    'Loop
    'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
    Pict = Form1.Picture7.hDC

    If GoOn = 1 Then
        FloodFill Pict, dx, dy, hBlue + 1
    End If
    
    Form1.Picture7.DrawMode = 13
    Form1.Picture7.FillStyle = 1
    
    
    
Next z



If z = 0 Then
    Form1.Picture7.ForeColor = red
    SP = 10
    EP = 13
ElseIf z = 1 Then
    Form1.Picture7.ForeColor = green
    SP = SP + 3
    EP = EP + 3
Else
    SP = SP + 3
    EP = EP + 3
    Form1.Picture7.ForeColor = blue
End If
Pict = Form1.Picture7.hDC

Form1.Picture7.ForeColor = green
SP = 10
EP = 13
For X = 0 To LenXoverSeqS1
    MoveToEx Pict, (30 + XDiffposS1(X) * XFactor), SP, PntAPI
    LineTo Pict, (30 + XDiffposS1(X) * XFactor), EP
Next X

Form1.Picture7.ForeColor = blue
SP = SP + 3
EP = EP + 3
For X = 0 To LenXoverSeqS2
    MoveToEx Pict, (30 + XDiffposS2(X) * XFactor), SP, PntAPI
    LineTo Pict, (30 + XDiffposS2(X) * XFactor), EP
Next X

Form1.Picture7.ForeColor = red
SP = SP + 3
EP = EP + 3
For X = 0 To LenXoverSeqS3
    MoveToEx Pict, (30 + XDiffposS3(X) * XFactor), SP, PntAPI
    LineTo Pict, (30 + XDiffposS3(X) * XFactor), EP
Next X
Form1.Picture7.DrawWidth = 2
Form1.Picture7.ForeColor = green
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + XDiffposS1(0) * XFactor, (PicHeight - (15 + ((XoverSeqNumTS1(0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
For X = 1 To LenXoverSeqS1
    LineTo Pict, 30 + XDiffposS1(X) * XFactor + XFactor, PicHeight - (15 + (((XoverSeqNumTS1(X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X

Form1.Picture7.ForeColor = blue
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + XDiffposS2(0) * XFactor, (PicHeight - (15 + ((XoverSeqNumTS2(0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
For X = 1 To LenXoverSeqS2
    LineTo Pict, 30 + XDiffposS2(X) * XFactor + XFactor, PicHeight - (15 + (((XoverSeqNumTS2(X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X

Form1.Picture7.ForeColor = red
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + XDiffposS3(0) * XFactor, (PicHeight - (15 + ((XoverSeqNumTS3(0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI
For X = 1 To LenXoverSeqS3
    LineTo Pict, 30 + XDiffposS3(X) * XFactor + XFactor, PicHeight - (15 + (((XoverSeqNumTS3(X) - Min) / (Max - Min))) * (PicHeight - 35))
Next 'X
Form1.Picture7.DrawWidth = 1

'Write sequenc names
Call WriteNames2(tseq1, tseq2, tseq3, red, green, blue)


 


Seq1 = tseq1
Seq2 = tseq2
Seq3 = tseq3

 
 

 ProbY = 1
 
 
 

 If BE < EN Then
     ProbX = (BE + (EN - BE) / 2)
 Else
     If EN > Len(StrainSeq(0)) - BE Then
         ProbX = EN / 2
     Else
         ProbX = BE + (Len(StrainSeq(0)) - BE) / 2
     End If
     
 End If
             
 
 EN = SuperEventList(XOverList(RelX, RelY).Eventnumber)
 ProbTest = dPValue
 Call PrintProbability
 MinPA = dPValue
 If Confirm(EN, 0) > 0 Then
    PT = 10 ^ (-ConfirmP(EN, 0))
Else
    PT = 1
End If
If MinPA < 0.9 Or MinPA <= LowestProb Then
     If (Confirm(EN, 8) = 1 Or Confirm(EN, 8) = 0) Or (XOverList(RelX, RelY).ProgramFlag <> 8 And XOverList(RelX, RelY).ProgramFlag <> 8 + AddNum) Then
    
            If ((Confirm(EN, 8) = 0 Or (Confirm(EN, 8) = 1) And MinPA < PT)) And MinPA < 1 And MinPA > 0 Then
                Confirm(EN, 8) = 1
                ConfirmP(EN, 8) = -Log10(MinPA)
                DoEvents
    
                If Form1.HScroll3.Value = 0 Then
                    Form1.HScroll3.Value = 1
                Else
                    Form1.HScroll3.Value = 0
                End If
    
            End If
    
        End If
 End If
 
 
 
 
 

    
    
    
If LenXoverSeqS1 >= LenXoverSeqS2 And LenXoverSeqS1 >= LenXoverSeqS3 Then
    LenXoverSeq = LenXoverSeqS1
ElseIf LenXoverSeqS2 >= LenXoverSeqS1 And LenXoverSeqS2 >= LenXoverSeqS3 Then
    LenXoverSeq = LenXoverSeqS2
ElseIf LenXoverSeqS3 >= LenXoverSeqS2 And LenXoverSeqS3 >= LenXoverSeqS1 Then
    LenXoverSeq = LenXoverSeqS3
End If
    
    
    'get everything into standard save/copy format
GPrintNum = 5 'six lines
NSites = LenXoverSeq + 1
ReDim GPrint(GPrintNum, NSites * 2 + 2), GPrintCol(GPrintNum), GPrintPos(GPrintNum, NSites * 2 + 2)

ReDim GVarPos(2, NSites)
ReDim Preserve PermRangeS1(1, LenXoverSeq + 1), PermRangeS2(1, LenXoverSeq + 1), PermRangeS3(1, LenXoverSeq + 1)
For X = LenXoverSeqS1 To NSites
    XDiffposS1(X) = XDiffposS1(LenXoverSeqS1)
    XoverSeqNumTS1(X) = XoverSeqNumTS1(LenXoverSeqS1)
    PermRangeS1(0, X) = PermRangeS1(0, LenXoverSeqS1)
    PermRangeS1(1, X) = PermRangeS1(1, LenXoverSeqS1)
Next X
For X = LenXoverSeqS2 To NSites
    XDiffposS2(X) = XDiffposS2(LenXoverSeqS2)
    XoverSeqNumTS2(X) = XoverSeqNumTS2(LenXoverSeqS2)
    PermRangeS2(0, X) = PermRangeS2(0, LenXoverSeqS2)
    PermRangeS2(1, X) = PermRangeS2(1, LenXoverSeqS2)
Next X
For X = LenXoverSeqS3 To NSites
    XDiffposS3(X) = XDiffposS3(LenXoverSeqS3)
    XoverSeqNumTS3(X) = XoverSeqNumTS3(LenXoverSeqS3)
    PermRangeS3(0, X) = PermRangeS3(0, LenXoverSeqS3)
    PermRangeS3(1, X) = PermRangeS3(1, LenXoverSeqS3)
Next X

For X = 1 To NSites
    GVarPos(0, X) = XDiffposS1(X - 1)
    GVarPos(1, X) = XDiffposS2(X - 1)
    GVarPos(2, X) = XDiffposS3(X - 1)
Next X


ReDim GCritval(0)


GLegend = "Height"
GPrintLen = NSites * 2 + 2 'how many points to plot
GPrintCol(0) = green 'line is black
GPrintCol(3) = hGreen 'line is grey

GPrintCol(1) = blue 'line is black
GPrintCol(4) = hBlue 'line is grey

GPrintCol(2) = red 'line is black
GPrintCol(5) = hRed 'line is grey

GPrintType = 0 'a normal line plot
GPrintMin(0) = Min   'bottom val
GPrintMin(1) = Max  'upper val



For X = 0 To NSites - 1
    
    GPrint(0, X) = XoverSeqNumTS1(X) 'GraphPlt(0, X)
    GPrint(0, NSites * 2 - X) = XoverSeqNumTS1(X)  'GraphPlt(0, X)
    GPrintPos(0, X) = XDiffposS1(X) 'PltPos(X)
    GPrintPos(0, NSites * 2 - X) = XDiffposS1(X) ' PltPos(X)
    
    GPrint(1, X) = XoverSeqNumTS2(X) 'GraphPlt(0, X)
    GPrint(1, NSites * 2 - X) = XoverSeqNumTS2(X)  'GraphPlt(0, X)
    GPrintPos(1, X) = XDiffposS2(X) 'PltPos(X)
    GPrintPos(1, NSites * 2 - X) = XDiffposS2(X) ' PltPos(X)
    
    GPrint(2, X) = XoverSeqNumTS3(X) 'GraphPlt(0, X)
    GPrint(2, NSites * 2 - X) = XoverSeqNumTS3(X)  'GraphPlt(0, X)
    GPrintPos(2, X) = XDiffposS3(X) 'PltPos(X)
    GPrintPos(2, NSites * 2 - X) = XDiffposS3(X) ' PltPos(X)
    
Next X


For X = 0 To NSites - 1
    'PValMap(DN, PermutationX * 0.01)
    GPrint(3, X) = PermRangeS1(0, X) 'GraphPlt(1, X)
    GPrint(3, NSites * 2 - X) = PermRangeS1(1, X) 'GraphPlt(2, X)
    GPrintPos(3, X) = XDiffposS1(X) 'PltPos(X)
    GPrintPos(3, NSites * 2 - X) = XDiffposS1(X) 'PltPos(X)
    
    GPrint(4, X) = PermRangeS2(0, X) 'GraphPlt(1, X)
    GPrint(4, NSites * 2 - X) = PermRangeS2(1, X) 'GraphPlt(2, X)
    GPrintPos(4, X) = XDiffposS2(X) 'PltPos(X)
    GPrintPos(4, NSites * 2 - X) = XDiffposS2(X) 'PltPos(X)
    
    GPrint(5, X) = PermRangeS3(0, X) 'GraphPlt(1, X)
    GPrint(5, NSites * 2 - X) = PermRangeS3(1, X) 'GraphPlt(2, X)
    GPrintPos(5, X) = XDiffposS3(X) 'PltPos(X)
    GPrintPos(5, NSites * 2 - X) = XDiffposS3(X) 'PltPos(X)
    
Next X
GPrintPos(3, GPrintLen - 1) = GPrintPos(3, 0)
GPrintPos(4, GPrintLen - 1) = GPrintPos(4, 0)
GPrintPos(5, GPrintLen - 1) = GPrintPos(5, 0)

GPrintPos(0, GPrintLen - 1) = GPrintPos(0, 0)
GPrintPos(0, GPrintLen) = GPrintPos(0, 0)
GPrint(0, GPrintLen - 1) = GPrint(0, 0)

GPrintPos(1, GPrintLen - 1) = GPrintPos(1, 0)
GPrintPos(1, GPrintLen) = GPrintPos(1, 0)
GPrint(1, GPrintLen - 1) = GPrint(1, 0)

GPrintPos(2, GPrintLen - 1) = GPrintPos(2, 0)
GPrintPos(2, GPrintLen) = GPrintPos(2, 0)
GPrint(2, GPrintLen - 1) = GPrint(2, 0)


GPrintPos(3, GPrintLen) = GPrintPos(3, 0)
GPrintPos(4, GPrintLen) = GPrintPos(4, 0)
GPrintPos(5, GPrintLen) = GPrintPos(5, 0)
GPrint(3, GPrintLen - 1) = GPrint(3, GPrintLen - 2)
GPrint(3, GPrintLen) = GPrint(3, 0)
GPrint(4, GPrintLen - 1) = GPrint(4, GPrintLen - 2)
GPrint(4, GPrintLen) = GPrint(4, 0)
GPrint(5, GPrintLen - 1) = GPrint(5, GPrintLen - 2)
GPrint(5, GPrintLen) = GPrint(5, 0)

GExtraTNum = 5
ReDim GExtraText(5)
GExtraText(0) = "Hyper-geometric random walk (" + StraiName(Seq1) + ")"
GExtraText(1) = "Hyper-geometric random walk (" + StraiName(Seq2) + ")"
GExtraText(2) = "Hyper-geometric random walk (" + StraiName(Seq3) + ")"

GExtraText(3) = "Bounds of 100 permuted datasets for " + StraiName(Seq1)
GExtraText(4) = "Bounds of 100 permuted datasets for " + StraiName(Seq2)
GExtraText(5) = "Bounds of 100 permuted datasets for " + StraiName(Seq3)

  Seq1 = tseq1
 Seq2 = tseq2
 Seq3 = tseq3
    
    
End Sub
Public Sub TSXOverB()
Dim BE As Long, EN As Long, BE2 As Long, EN2 As Long, nK As Long, nM As Long, nN As Long, nL As Long, X As Long, CurrentHeight As Long, XoverSeqNumTS() As Long, XOverSeqNumPerms() As Long, Y As Long, MaxSeen As Long, MaxDescentSeen As Long
Dim ii As Long, dPValue As Double


tseq1 = Seq1
tseq2 = Seq2
tseq3 = Seq3
XX = RelX
XX = RelY

If XTableFlag = 0 Then Call Build3SeqTable
If XOverList(RelX, RelY).ProgramFlag = 8 Or XOverList(RelX, RelY).ProgramFlag = 8 + AddNum Then
    Seq3 = CLng(Abs(XOverList(RelX, RelY).DHolder))
    'XOverList(RelX, RelY).DHolder = 51
    'Seq1 = 6
    'XOverList(RelX, RelY).LHolder = 75
    If Seq3 <> tseq1 And Seq3 <> tseq2 And Seq3 <> tseq3 Then
        For X = 0 To Nextno
            If TreeTrace(X) = Seq3 And (X = tseq1 Or X = tseq2 Or X = tseq3) Then
                Seq3 = X: Exit For
            End If
        Next X
    End If
    If Seq3 = tseq1 Then
        Seq2 = tseq2
        Seq1 = tseq3
    ElseIf Seq3 = tseq2 Then
        Seq2 = tseq1
        Seq1 = tseq3
    ElseIf Seq3 = tseq3 Then
        Seq2 = tseq1
        Seq1 = tseq2
    Else
        
        
        Seq1 = tseq1
        Seq2 = tseq2
        Seq3 = tseq3
        Call TSXOverC
        Exit Sub
    End If
Else
    Seq1 = tseq3
    Seq2 = tseq1
    Seq3 = tseq2
End If
XX = StraiName(Seq3)



'Call TSXOver(0)
ReDim XoverSeqNumTS(Len(StrainSeq(0)))
ReDim XOverSeqNumPerms(Len(StrainSeq(0)), 100)

Dim tXOverSeqNum() As Long
ReDim tXOverSeqNum(Len(StrainSeq(0)))
Rnd (GetTickCount)
Dim Tmp As Long, THold(2) As Long, SeqRnd() As Integer, dPValue2 As Double, xPosdiffx() As Long, xDiffposx() As Long, sn As Long, sm As Long, sk As Long, sbe As Long, sen As Long
XX = UBound(SeqNum, 2)
ReDim SeqRnd(Len(StrainSeq(0)), 2)
        
        BE = 0
        EN = 0
        nM = 0
        nN = 0
        nK = 0
        LenXoverSeq = FindSubSeqTS2(Len(StrainSeq(0)), Seq1, Seq2, Seq3, BE, EN, BE2, EN2, nM, nN, nK, nL, XPosDiff(0), XDiffpos(0), SeqNum(0, 0), XoverSeqNumTS(0), MissingData(0, 0))
        
        '2,7,3:2642,1208,1208,2642:26,20,13,19
        '1312,2309
        '2352,1306
        
        '1306,2309
        '2309,1306
        
    
        
        Call CheckWrap(LenXoverSeq, nK, BE, EN, 1, XDiffpos(), XPosDiff(), XoverSeqNumTS())
        Call CheckWrap(LenXoverSeq, nL, BE2, EN2, -1, XDiffpos(), XPosDiff(), XoverSeqNumTS())
        '13,2653,1208
        '19,1214,2642
       
        GetTSPVal WF, nM, nN, nK, dPValue
        GetTSPVal WF2, nN, nM, nL, dPValue2
        MultNeg = 0
        If dPValue2 < dPValue Then
            Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
        End If
        
        If MultNeg = 1 Then
            For X = 0 To LenXoverSeq
                XoverSeqNumTS(X) = XoverSeqNumTS(X) * -1
            Next X
        End If
            '1.646-4,699,18
            
        Call CheckSplit3Seq(BE, EN, nM, nN, dPValue, LenXoverSeq, XoverSeqNumTS(), XPosDiff(), XDiffpos())
        If dPValue > dPValue2 Or FindallFlag = 1 Then
            If MultNeg = 0 Then
                For X = 0 To LenXoverSeq
                    XoverSeqNumTS(X) = XoverSeqNumTS(X) * -1
                Next X
                'Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
            End If
            Call CheckSplit3Seq(BE2, EN2, nN, nM, dPValue2, LenXoverSeq, XoverSeqNumTS(), XPosDiff(), XDiffpos())
            If FindallFlag = 0 And dPValue > dPValue2 Then
                If MultNeg = 1 Then
                    MultNeg = 0
                    'For X = 0 To LenXoverSeq
                    '    XoverSeqNumTS(X) = XoverSeqNumTS(X) * -1
                    'Next X
                Else
                    MultNeg = 1
                    'For X = 0 To LenXoverSeq
                    '    XoverSeqNumTS(X) = XoverSeqNumTS(X) * -1
                    'Next X
                    
                End If
                Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
                
            End If
        
        End If
        
        
        
        
        If BE = 0 Then BE = XDiffpos(0)
        If LenXoverSeq < 1 Then Exit Sub
        
        
        Dim tMissingData() As Long
        ReDim tMissingData(Len(StrainSeq(0)), 2)
        For X = 0 To Len(StrainSeq(0))
            tMissingData(X, 0) = MissingData(X, Seq1)
            tMissingData(X, 1) = MissingData(X, Seq2)
            tMissingData(X, 2) = MissingData(X, Seq3)
        Next X
        
        xtSeq1 = Seq1
        xtSeq2 = Seq2
        
        If MultNeg = 1 Then
            Seq1 = xtSeq2
            Seq2 = xtSeq1
        End If
        
        For X = 1 To 100
            ReDim xPosdiffx(Len(StrainSeq(0))), xDiffposx(Len(StrainSeq(0)))
            For z = 1 To Len(StrainSeq(0))
                
                SeqRnd(z, 0) = SeqNum(z, Seq1)
                SeqRnd(z, 1) = SeqNum(z, Seq2)
                SeqRnd(z, 2) = SeqNum(z, Seq3)
            Next z
           
            For z = 1 To Len(StrainSeq(0))
                NewPos = CLng(Rnd * (Len(StrainSeq(0)) - 1)) + 1
                If tMissingData(z, 0) = 0 And tMissingData(z, 1) = 0 And tMissingData(z, 2) = 0 And tMissingData(NewPos, 0) = 0 And tMissingData(NewPos, 1) = 0 And tMissingData(NewPos, 2) = 0 Then
                    THold(0) = SeqRnd(z, 0)
                    THold(1) = SeqRnd(z, 1)
                    THold(2) = SeqRnd(z, 2)
                    SeqRnd(z, 0) = SeqRnd(NewPos, 0)
                    SeqRnd(z, 1) = SeqRnd(NewPos, 1)
                    SeqRnd(z, 2) = SeqRnd(NewPos, 2)
                    SeqRnd(NewPos, 0) = THold(0)
                    SeqRnd(NewPos, 1) = THold(1)
                    SeqRnd(NewPos, 2) = THold(2)
                End If
            Next z
            ReDim tXOverSeqNum(Len(StrainSeq(0)))
            sbe = 0
            sen = 0
            sm = 0
            sn = 0
            sk = 0
            
            LenXoverSeqx = FindSubSeqTS2(Len(StrainSeq(0)), 0, 1, 2, sbe, sen, BE2, EN2, sm, sn, sk, nL, xPosdiffx(0), xDiffposx(0), SeqRnd(0, 0), tXOverSeqNum(0), tMissingData(0, 0))
            For z = 0 To Len(StrainSeq(0))
                XOverSeqNumPerms(z, X) = tXOverSeqNum(z)
            Next z
        Next X
        
        Seq1 = xtSeq1
        Seq2 = xtSeq2
   
    
    'If nN > 0 And nK = 1 Then Exit Sub
    
    'If nN - nM = nK Then Exit Sub
    
    'Compute P-vale for maxdescent=(nK)
    
   ' if( m < mSize && n < nSize && k < kSize && n < kSize && n < jSize )
   ' {
   '     double dPValue = ((double)0);
    
    
    
    
    
    
    
   ' performs a Dunn-Sidak correction for pval with m trials
    xpvalue = dPValue * MCCorrection
    If dPValue >= 1 Then
        dPValue = 1
    Else
        dPValue = 1 - (1 - dPValue) ^ MCCorrection
    End If
    'characterise the event
    If dPValue = 0 Then dPValue = xpvalue
    Dim Max As Double, Min As Double
    Max = -LenXoverSeq
    Min = LenXoverSeq
    For X = 0 To LenXoverSeq
        If Max < XoverSeqNumTS(X) Then Max = XoverSeqNumTS(X)
        If Min > XoverSeqNumTS(X) Then Min = XoverSeqNumTS(X)
        
    
    Next X
    
    For z = 1 To 100
        For X = 0 To LenXoverSeq
            If Max < XOverSeqNumPerms(X, z) Then Max = XOverSeqNumPerms(X, z)
            If Min > XOverSeqNumPerms(X, z) Then Min = XOverSeqNumPerms(X, z)
            
        
        Next X
    Next z
    
    Dim YScaleFactor As Double, PntAPI As POINTAPI, Pict As Long
    Form1.Picture7.Picture = LoadPicture()
    YScaleFactor = 0.85
    PicHeight = Form1.Picture7.Height * YScaleFactor
    
    XFactor = ((Form1.Picture7.Width - 40) / Len(StrainSeq(0)))
    
    Call DoAxes(0, Len(StrainSeq(0)), -1, Max, Min, 1, "Height")
    Pict = Form1.Picture7.hDC
    
    Form1.Picture7.ForeColor = RGB(127, 127, 127)
    Dim PermRange() As Long
    ReDim PermRange(1, LenXoverSeq)
    
    
    For X = 0 To LenXoverSeq
        PermRange(1, X) = LenXoverSeq
        PermRange(0, X) = -LenXoverSeq
        For Y = 1 To 100
            If XOverSeqNumPerms(X, Y) > PermRange(0, X) Then PermRange(0, X) = XOverSeqNumPerms(X, Y)
            If XOverSeqNumPerms(X, Y) < PermRange(1, X) Then PermRange(1, X) = XOverSeqNumPerms(X, Y)
        Next Y
    Next X
    
    MoveToEx Pict, 30 + XDiffpos(0) * XFactor, (PicHeight - (15 + ((PermRange(0, 0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI

    For X = 0 To LenXoverSeq
        
        LineTo Pict, 30 + XDiffpos(X) * XFactor + XFactor, PicHeight - (15 + (((PermRange(1, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    
    For X = LenXoverSeq To 0 Step -1
        
        LineTo Pict, 30 + XDiffpos(X) * XFactor + XFactor, PicHeight - (15 + (((PermRange(0, X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X
    LineTo Pict, 30 + XDiffpos(0) * XFactor + XFactor, PicHeight - (15 + (((PermRange(1, 0) - Min) / (Max - Min))) * (PicHeight - 35))
    
    'MoveToEx Pict, 30 + XDiffpos(2) * XFactor + XFactor, 0, PntAPI
    'LineTo Pict, 30 + XDiffpos(2) * XFactor + XFactor, 1000
    Form1.Picture7.FillStyle = 0
    Form1.Picture7.FillColor = RGB(128, 128, 128)
    Pict = Form1.Picture7.hDC
    X = 0
    GoOn = 0
    Dim MaxR As Long, WinR As Long
    WinR = -1
    
    MaxR = 0
    Do While X <= LenXoverSeq
        If PermRange(0, X) - PermRange(1, X) > MaxR Then
            MaxR = PermRange(0, X) - PermRange(1, X)
            WinR = X
        End If
        X = X + 1
    Loop
    
    X = WinR
    If X < LenXoverSeq Then
        dy = PicHeight - (15 + ((((PermRange(1, X) + (PermRange(0, X) - PermRange(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
        dx = 30 + XDiffpos(X) * XFactor + XFactor
    Else
        dy = PicHeight - (15 + ((((PermRange(1, X) + (PermRange(0, X) - PermRange(1, X)) / 2) - Min) / (Max - Min))) * (PicHeight - 35))
        dx = 30 + XDiffpos(X) * XFactor + XFactor
    End If
    GoOn = 1
       
   
    'ExtFloodFill Pict, XDiffpos(0) * XFactor + XFactor + 1, dy, RGB(128, 128, 128), 1
    Pict = Form1.Picture7.hDC

    If GoOn = 1 Then
        FloodFill Pict, dx, dy, RGB(127, 127, 127)
    End If
    
    'LineTo Pict, 30 + XDiffpos(0) * XFactor + XFactor + 2, dy
    Form1.Picture7.FillStyle = 1
    If TManFlag <> 22 Then
        Call Highlight
    Else
        BE = 0
        EN = Len(StrainSeq(0))
        GBlockNum = -1
    End If
    Form1.Picture7.ForeColor = 0
    
    
    
    MoveToEx Pict, 30 + XDiffpos(0) * XFactor, (PicHeight - (15 + ((XoverSeqNumTS(0) - Min) / (Max - Min)) * (PicHeight - 35))), PntAPI

    For X = 1 To LenXoverSeq
    
        LineTo Pict, 30 + XDiffpos(X) * XFactor + XFactor, PicHeight - (15 + (((XoverSeqNumTS(X) - Min) / (Max - Min))) * (PicHeight - 35))
    Next 'X

    Seq1 = tseq1
    Seq2 = tseq2
    Seq3 = tseq3

    
    ProbY = 1
    If BE < EN Then
        ProbX = (BE + (EN - BE) / 2)
    Else
        If EN > Len(StrainSeq(0)) - BE Then
            ProbX = EN / 2
        Else
            ProbX = BE + (Len(StrainSeq(0)) - BE) / 2
        End If
        
    End If
                
    
    
    ProbTest = dPValue
    Call PrintProbability
    
    If TManFlag <> 22 Then
        ReDim CXoverSeq(2)
        
        
        For X = 1 To LenXoverSeq
            CXoverSeq(0) = CXoverSeq(0) + Mid(StrainSeq(TreeTrace(Seq1)), XDiffpos(X), 1)
            CXoverSeq(1) = CXoverSeq(1) + Mid(StrainSeq(TreeTrace(Seq2)), XDiffpos(X), 1)
            CXoverSeq(2) = CXoverSeq(2) + Mid(StrainSeq(TreeTrace(Seq3)), XDiffpos(X), 1)
        Next X
    End If
    Seq1 = tseq1
    Seq2 = tseq2
    Seq3 = tseq3
    
    
    'get everything into standard save/copy format
GPrintNum = 1 'two lines
NSites = LenXoverSeq + 1
ReDim GPrint(GPrintNum, NSites * 2 + 2), GPrintCol(GPrintNum), GPrintPos(GPrintNum, NSites * 2 + 2)

ReDim GVarPos(0, NSites)
For X = 1 To NSites
    GVarPos(0, X) = XDiffpos(X - 1)
Next X


ReDim GCritval(0)


GLegend = "Height"
GPrintLen = NSites * 2 + 2 'how many points to plot
GPrintCol(0) = 0 'line is black
GPrintCol(1) = RGB(128, 128, 128) 'line is grey


GPrintType = 0 'a normal line plot
GPrintMin(0) = Min   'bottom val
GPrintMin(1) = Max  'upper val



For X = 0 To NSites - 1
    
    GPrint(0, X) = XoverSeqNumTS(X) 'GraphPlt(0, X)
    GPrint(0, NSites * 2 - X) = XoverSeqNumTS(X)  'GraphPlt(0, X)
    GPrintPos(0, X) = XDiffpos(X) 'PltPos(X)
    GPrintPos(0, NSites * 2 - X) = XDiffpos(X) ' PltPos(X)
    
    
Next X


For X = 0 To NSites - 1
    'PValMap(DN, PermutationX * 0.01)
    GPrint(1, X) = PermRange(0, X) 'GraphPlt(1, X)
    GPrint(1, NSites * 2 - X) = PermRange(1, X) 'GraphPlt(2, X)
    GPrintPos(1, X) = XDiffpos(X) 'PltPos(X)
    GPrintPos(1, NSites * 2 - X) = XDiffpos(X) 'PltPos(X)
    
    
Next X
GPrintPos(1, GPrintLen - 1) = GPrintPos(1, 0)

GPrintPos(0, GPrintLen - 1) = GPrintPos(0, 0)
GPrintPos(1, GPrintLen) = GPrintPos(1, 0)
GPrintPos(0, GPrintLen) = GPrintPos(0, 0)

GPrint(0, GPrintLen - 1) = GPrint(0, 0)
GPrint(1, GPrintLen - 1) = GPrint(1, GPrintLen - 2)
GPrint(1, GPrintLen) = GPrint(1, 0)

GExtraTNum = 1
ReDim GExtraText(1)
GExtraText(0) = "Hyper-geometric random walk"
GExtraText(1) = "Bounds of 100 permuted datasets"

    
    
    
End Sub
Public Sub TSChecking()
    Dim WeightMod() As Long, Scratch() As Integer
    'This is the subroutine that is called when checking RDP results
    Screen.MousePointer = 11
    OptFlag = -1
    If CurrentCheck = 0 Then
        xSpacerFlag = SpacerFlag
        XOverWindowX = CDbl(Form3.Text2.Text)
        
        Call XOverIII(0)
        
        SpacerFlag = xSpacerFlag
        Call FindSubSeqChi
    ElseIf CurrentCheck = 1 Then
        
        'If togglex = 0 Then
            Call GCCheck(0)
            Call FindSubSeqChi
        'Else
        '    Call GCXoverE
        'End If
        'togglex = togglex + 1
        'If togglex = 2 Then togglex = 0
        
    ElseIf CurrentCheck = 13 Then
        Call PXoverD
       Call FindSubSeqChi
     ElseIf CurrentCheck = 16 Then
        Call TSXOverC
        
    ElseIf CurrentCheck = 2 Then
        
        Call FindSubSeqBS
        
        If DoScans(0, 2) = 1 And UBound(BSFilePos, 2) > 0 Then

            Call BSXoverL(0)

        Else

            
                
                s1col = Yellow
                s1colb = LYellow
                s2col = Purple
                s2colb = LPurple
                s3col = green
                s2colb = LGreen
                Call FindSubSeqBS
                
                ReDim Scratch(BSStepWin), WeightMod(BSBootReps, BSStepWin - 1)
                Dummy = SEQBOOT2(BSRndNumSeed, BSBootReps, BSStepWin, Scratch(0), WeightMod(0, 0))
                Call BSXoverM(0, 0, WeightMod())

               

            

        End If
        
        Call FindSubSeqChi

    ElseIf CurrentCheck = 3 Then
        
        Call FindSubSeqBS
        
        
            s1col = Yellow
            s1colb = LYellow
            s2col = Purple
            s2colb = LPurple
            s3col = green
            s2colb = LGreen
            Call FindSubSeqBS
            
                ReDim Scratch(BSStepWin), WeightMod(BSBootReps, BSStepWin - 1)
                Dummy = SEQBOOT2(BSRndNumSeed, BSBootReps, BSStepWin, Scratch(0), WeightMod(0, 0))
                Call BSXoverM(0, 0, WeightMod())


        Call FindSubSeqChi
        
        
    ElseIf CurrentCheck = 4 Then
        
        Call MCXoverG(0)
        Call FindSubSeqChi
       
    ElseIf CurrentCheck = 5 Then
        
        Call SSXoverB(0)
       
       Call FindSubSeqChi
       
    ElseIf CurrentCheck = 6 Then

       
            
            Call FindSubSeqBS
            
            Call LXoverB(0, 1)
            Call FindSubSeqChi
          

    ElseIf CurrentCheck = 7 Then

  

    ElseIf CurrentCheck = 8 Then
        
        Call FindSubSeqBS
        
        Call DXoverE
        Call FindSubSeqChi
    ElseIf CurrentCheck = 9 Then

        
            
            Call FindSubSeqBS
            
            Call TXover3
            Call FindSubSeqChi

    ElseIf CurrentCheck = 10 Then
        
        Call CXoverC(0)
        Call FindSubSeqChi
    ElseIf CurrentCheck = 11 Then

        Call RecOverview2

    ElseIf CurrentCheck = 12 Then
        Call RecombMap
    ElseIf CurrentCheck = 14 Then
        Call RecombMapII
    ElseIf CurrentCheck = 15 Then
        Call RecombMapIII
    End If
    Dim AH(2) As Long

    Screen.MousePointer = 0
End Sub

Public Function Get3SeqPval(M, N, K, J)
'VB version of Boni's float ytable::prob(int m, int n, int k, int j) in ProbTables.ccp
If M > UBound(YTable, 1) - 1 Or N > UBound(YTable, 1) - 1 Then Get3SeqPval = 1: Exit Function
If J > UBound(YTable, 3) - 1 Or K > UBound(YTable, 3) - 1 Then Get3SeqPval = 1: Exit Function
Dim d As Single
    
    
    d = -1
    If J > K Then
    
        Get3SeqPval = 0
        Exit Function
    ElseIf K > N Or K < N - M Or J > N Or J < N - M Then
    
        Get3SeqPval = 0
        Exit Function
    ElseIf N = 0 Then
    
        If K = 0 And J = 0 Then
        
            Get3SeqPval = 1
            Exit Function
        
        Else
            Get3SeqPval = 0
            Exit Function
        End If
    ElseIf M = 0 Then
   
        If K = N And J = N Then
       
            Get3SeqPval = 1
            Exit Function
        
        Else
       
            Get3SeqPval = 0
            Exit Function
       End If
     ElseIf K = 0 And J = 0 Then
    
        If N = 0 Then
        
            Get3SeqPval = 1
            Exit Function
        
        Else
        
            Get3SeqPval = 0
            Exit Function
        End If
    End If
    
    
    
'    if( metatable[ m*nSize + n ] == NULL )
'    {
'        metatable[ m*nSize + n ] = new float[ jSize*kSize ];
'    T = metatable[ m*nSize + n ];
    
'    // now fill it up with (-1)s
'    for(a=0; a < kSize; a++)
'        for(b=0; b < jSize; b++ )
'        {
'            T[ jSize*a + b ] = ((float)-1);
'        }
'
'    // and make sure you record that you have allocated memory for
'    // another pointer in the metatable
'    nNumAllocatedInMeta++;
'    }
'    Else
'    {
'        T = metatable[ m*nSize + n ];
'    }
    
    ' if it's in the table, just return it
    
    
   
        
        
        If YTable(M, N, K, J) >= 0 Then
            
            Get3SeqPval = YTable(M, N, K, J)
            '50,1,1,1
            '50,1,1,0
        Else
            If J = 0 Then
            
                d = (M / (N + M)) * (Get3SeqPval(M - 1, N, K, 1) + Get3SeqPval(M - 1, N, K, 0))
                
            Else 'j > 0 Then
                If K = J Then
                    d = (N / (N + M)) * (Get3SeqPval(M, N - 1, J - 1, J - 1) + Get3SeqPval(M, N - 1, J, J - 1))
                Else ' k > j Then
                    d = (M / (N + M)) * Get3SeqPval(M - 1, N, K, J + 1) + (N / (N + M)) * Get3SeqPval(M, N - 1, K, J - 1)
                End If
            End If
            
            
            YTable(M, N, K, J) = d
            
            Get3SeqPval = d
        End If
    
    
End Function

Public Sub SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
Dim T As Double
T = BE: BE = BE2: BE2 = T
T = EN: EN = EN2: EN2 = T

Tmp = nM
nM = nN
nN = Tmp

T = nK: nK = nL: nL = T
T = dPValue: dPValue = dPValue2: dPValue2 = T
MultNeg = 1
End Sub
Public Sub CheckSplit3Seq(BE, EN, nM, nN, dPValue, LenXoverSeq, XoverSeqNumTS() As Long, XPosDiff() As Long, XDiffpos() As Long)
Dim SplitFlag As Byte, tdPValue As Double, tBE As Long, tEN As Long


tBE = BE
tEN = EN
X = BE - 1
Do While X <> EN
    If X > Len(StrainSeq(0)) Then
        X = 1
        If EN = 0 Or EN = 1 Then
            Exit Do
        End If
    End If
    If X < 0 Then X = Len(StrainSeq(0))
    If MissingData(X, Seq1) = 1 Or MissingData(X, Seq2) = 1 Or MissingData(X, Seq3) = 1 Then
        '8.799e-5'26,47,33
        
        dPValue = 1
        tdPValue = SubPVal(nM, nN, LenXoverSeq, BE, X, XPosDiff(), XDiffpos, XoverSeqNumTS())
        If tdPValue < dPValue Then
            dPValue = tdPValue
            tBE = BE
            tEN = X - 1
        End If
        Exit Do
    End If
    If X = EN Then Exit Do
    X = X + 1
    If X = EN Then Exit Do
Loop

If X <> EN Then
    X = EN + 1
    Do While X <> BE
        If X < 0 Then X = Len(StrainSeq(0))
        If X > Len(StrainSeq(0)) Then X = 1
        If MissingData(X, Seq1) = 1 Or MissingData(X, Seq2) = 1 Or MissingData(X, Seq3) = 1 Then
            
            tdPValue = SubPVal(nM, nN, LenXoverSeq, X, EN, XPosDiff(), XDiffpos, XoverSeqNumTS())
            If tdPValue < dPValue Then
                dPValue = tdPValue
                tBE = X + 1
                tEN = EN
            End If
            Exit Do
        End If
        If X = BE Then Exit Do
        X = X - 1
    Loop

End If
EN = tEN
BE = tBE
End Sub

Public Sub CheckWrap(LenXoverSeq, nK, BE, EN, NegMod, XDiffpos() As Long, XPosDiff() As Long, XoverSeqNumTS() As Long)
    If BE = 0 Then BE = XDiffpos(0)
    
    If CircularFlag = 1 Or X = X Then
        'need to carry on and check to see if mindescent gets lower with wrapping
       MaxSeen = XoverSeqNumTS(XPosDiff(BE)) * NegMod
       tempAscent = XoverSeqNumTS(LenXoverSeq) * NegMod
       If BE < EN Then
            For X = 0 To XPosDiff(BE)
                 
                 If tempAscent + XoverSeqNumTS(X) * NegMod > MaxSeen Then
                    MaxSeen = tempAscent + XoverSeqNumTS(X) * NegMod
                    BE = XDiffpos(X)
                 End If
                 If MaxSeen - (tempAscent + XoverSeqNumTS(X) * NegMod) > nK Then
                    nK = MaxSeen - (tempAscent + XoverSeqNumTS(X) * NegMod)
                    EN = XDiffpos(X)
                End If
                 
             Next X
             
             'XX = MissingData(EN, Seq3)
        Else
            '0 to 15
            For X = 0 To XPosDiff(EN)
                 
                 '@'@'@'@'@'@'@'@
                 If tempAscent + XoverSeqNumTS(X) * NegMod > MaxSeen Then
                    MaxSeen = tempAscent + XoverSeqNumTS(X) * NegMod
                    BE = XDiffpos(X)
                 End If
                 If MaxSeen - (tempAscent + XoverSeqNumTS(X) * NegMod) > nK Then
                    nK = MaxSeen - (tempAscent + XoverSeqNumTS(X) * NegMod)
                    EN = XDiffpos(X)
                    X = X
                End If
                 
             Next X
        X = X
        End If
    
        
        
    End If
 '   If NegMod = 1 Then
        'If BE < EN Then 'ie highest point first, lowest point last
        
            If XPosDiff(BE) < LenXoverSeq Then
                BE = XDiffpos(XPosDiff(BE) + 1)
            Else
                BE = XDiffpos(0)
            End If
        'Else 'ie lowest point first, higest point last
        '    If XPosDiff(EN) < LenXOverSeq Then
        '        EN = XDiffPos(XPosDiff(EN) + 1)
        '    Else
        '        EN = XDiffPos(0)
        '    End If
        
        'End If
 '   Else
        'If BE > EN Then 'ie lowest point first, highest point last
    
 '           If XPosDiff(BE) < LenXOverSeq Then
 '               BE = XDiffPos(XPosDiff(BE) + 1)
 '           Else
 '               BE = XDiffPos(0)
 '           End If
        'Else 'ie highest point first, lowest point last
        '    If XPosDiff(EN) < LenXOverSeq Then
        '        EN = XDiffPos(XPosDiff(EN) + 1)
        '    Else
        '        EN = XDiffPos(0)
        '    End If
       '
       ' End If
 '   End If
    
    '@
    If CircularFlag = 0 Then
       
        If BE > EN Then
            If XPosDiff(EN) < LenXoverSeq Then
                TE = XDiffpos(XPosDiff(EN) + 1)
            Else
                TE = 1
            End If
            If XPosDiff(BE) > 0 Then
                EN = XDiffpos(XPosDiff(BE) - 1)
            Else
                EN = Len(StrainSeq(0))
            End If
             
            BE = TE
        End If
    
    End If
    
    
    
    Exit Sub
    'Centre the breakpoints
    If XPosDiff(EN) < LenXoverSeq Then
        EN = XDiffpos(XPosDiff(EN)) + (XDiffpos(XPosDiff(EN) + 1) - XDiffpos(XPosDiff(EN))) / 2
    Else
        EN = Len(StrainSeq(0))
    End If
    If XPosDiff(BE) > 0 Then
        BE = XDiffpos(XPosDiff(BE) - 1) + (XDiffpos(XPosDiff(BE)) - XDiffpos(XPosDiff(BE) - 1)) / 2
    Else
        BE = 1
    End If
End Sub
Public Sub TSXOver(FindallFlag)


Dim tempAscent As Long, BE As Long, EN As Long, BE2 As Long, EN2 As Long, nK As Long, nM As Long, nN As Long, nL As Long, X As Long, CurrentHeight As Long, Y As Long, MaxSeen As Long, MaxDescentSeen As Long
'3.610
    If XTableFlag = 0 Then Call Build3SeqTable
    '@'@'@'@'@'@'@'@'@'@'@
    
    
    If SEventNumber = 0 Then
        LenXoverSeq = FindSubSeqTS(Len(StrainSeq(0)), Seq1, Seq2, Seq3, BE, EN, BE2, EN2, nM, nN, nK, nL, XPosDiff(0), XDiffpos(0), SeqNum(0, 0), XoverSeqNumTS(0))
    Else
        LenXoverSeq = FindSubSeqTS2(Len(StrainSeq(0)), Seq1, Seq2, Seq3, BE, EN, BE2, EN2, nM, nN, nK, nL, XPosDiff(0), XDiffpos(0), SeqNum(0, 0), XoverSeqNumTS(0), MissingData(0, 0))
    End If
    If LenXoverSeq < 3 Then Exit Sub
    
'50.016
'43.235

' scan through the loop 5.391
'24.469 just to do currntheight
'34.704 - adding vals to xoverseqnum
    
    
    If X = X Then
        Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeq, 1, CircularFlag, nK, BE, EN, XDiffpos(0), XPosDiff(0), XoverSeqNumTS(0))
        Dummy = CheckwrapC(Len(StrainSeq(0)), LenXoverSeq, -1, CircularFlag, nL, BE2, EN2, XDiffpos(0), XPosDiff(0), XoverSeqNumTS(0))
    Else
        Call CheckWrap(LenXoverSeq, nK, BE, EN, 1, XDiffpos(), XPosDiff(), XoverSeqNumTS())
        Call CheckWrap(LenXoverSeq, nL, BE2, EN2, -1, XDiffpos(), XPosDiff(), XoverSeqNumTS())
    End If
    '26, 543-2309
    '11, 2352,398 :1775,398

    Dim ii As Long, dPValue As Double, MultNeg As Byte, T As Double
    
    'Calculate/read P-value for maxdescent=(nK)
    
    GetTSPVal WF, nM, nN, nK, dPValue
    GetTSPVal WF2, nN, nM, nL, dPValue2
    'dPValue2 = 10
    MultNeg = 0
    '@
    If dPValue2 < dPValue Then
        Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, dPValue, dPValue2, MultNeg)
    End If
    
    If nN > 0 And nK = 1 Then Exit Sub
    If nN - nM = nK Then Exit Sub
    
    Dim odPValue As Double, dPValueX As Double
    odPValue = dPValue
   ' performs a Dunn-Sidak correction for pval with m trials
    If MCFlag = 0 Then
        xpvalue = dPValue * MCCorrection
        If dPValue >= 1 Then
            dPValue = 1
        Else
            dPValue = 1 - (1 - dPValue) ^ MCCorrection
        End If
    Else
        xpvalue = dPValue
    End If
    
    'characterise the event
    If ((dPValue < 1 And dPValue <= LowestProb) Or (dPValue = 1 And xpvalue < LowestProb)) And xpvalue > 0 Then
    
        
        If SEventNumber > 0 Then
            If MultNeg = 1 Then
                For X = 0 To LenXoverSeq
                    XoverSeqNumTS(X) = XoverSeqNumTS(X) * -1
                Next X
            End If
            '1.646-4,699,18
            
            Call CheckSplit3Seq(BE, EN, nM, nN, odPValue, LenXoverSeq, XoverSeqNumTS(), XPosDiff(), XDiffpos())
            If odPValue > dPValue2 Or FindallFlag = 1 Then
                For X = 0 To LenXoverSeq
                    XoverSeqNumTS(X) = XoverSeqNumTS(X) * -1
                Next X
                Call CheckSplit3Seq(BE2, EN2, nN, nM, dPValue2, LenXoverSeq, XoverSeqNumTS(), XPosDiff(), XDiffpos())
                If FindallFlag = 0 And odPValue > dPValue2 Then
                    
                    Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, odPValue, dPValue2, MultNeg)
                    
                End If
            
            End If
            
            If MCFlag = 0 Then
                xpvalue = odPValue * MCCorrection
                
                If odPValue >= 1 Then
                    dPValue = 1
                Else
                    dPValue = 1 - (1 - odPValue) ^ MCCorrection
                End If
                
            Else
                xpvalue = odPValue
            End If
        End If
        If dPValue = 0 Then dPValue = xpvalue
        For z = 0 To 1
            If ((dPValue < 1 And dPValue <= LowestProb) Or (dPValue = 1 And xpvalue < LowestProb)) And xpvalue > 0 Then
                'Keep track of signal numbers
                If LongWindedFlag = 1 Then
                    If CurrentXOver(Seq1) <= CurrentXOver(Seq2) And CurrentXOver(Seq1) <= CurrentXOver(Seq3) Then
                        ActiveSeq = Seq1
                        ActiveMinorP = Seq2
                        ActiveMajorP = Seq3
                    ElseIf CurrentXOver(Seq2) <= CurrentXOver(Seq1) And CurrentXOver(Seq2) <= CurrentXOver(Seq3) Then
                        ActiveSeq = Seq2
                        ActiveMinorP = Seq1
                        ActiveMajorP = Seq3
                    Else
                        ActiveSeq = Seq3
                        ActiveMinorP = Seq1
                        ActiveMajorP = Seq2
                    End If
                End If
                oRecombNo(100) = oRecombNo(100) + 1
                oRecombNo(8) = oRecombNo(8) + 1
                
                If APermFlag = 0 Then
                    Call UpdateXOList3(ActiveSeq, CurrentXOver(), XOverList(), 8, dPValue, SIP)
                Else
                    SIP = 1
                End If
                'Exit Sub
                If SIP = -1 Then
                    If DoneRedo = 0 Then
                        DoneRedo = 1
                        Call AddToRedoList(8, Seq1, Seq2, Seq3)
                    End If
                Else
                    XOverList(ActiveSeq, SIP).Daughter = ActiveSeq
                    XOverList(ActiveSeq, SIP).MajorP = ActiveMajorP
                    XOverList(ActiveSeq, SIP).MinorP = ActiveMinorP
                    XOverList(ActiveSeq, SIP).SBPFlag = 0
                    XOverList(ActiveSeq, SIP).ProgramFlag = 8
                    BWarn = 0
                    EWarn = 0
                    Call CentreBP(BE, EN, XPosDiff(), XDiffpos(), BWarn, EWarn, 10, LenXoverSeq)
                    
                    XOverList(ActiveSeq, SIP).Beginning = BE
                    XOverList(ActiveSeq, SIP).Ending = EN
                    XOverList(ActiveSeq, SIP).Probability = dPValue
                    XOverList(ActiveSeq, SIP).PermPVal = WF
                    
                    XOverList(ActiveSeq, SIP).BeginP = 0
                    XOverList(ActiveSeq, SIP).EndP = 0
                    XOverList(ActiveSeq, SIP).DHolder = Seq3
                    
                    If ShortOutFlag = 1 Then
                        ShortOutput(8) = 1
                        AbortFlag = 1
                        Exit Sub
                    End If
                    
                    If XOverList(ActiveSeq, SIP).DHolder = 246 Then
                        '238,216,246
                        XX = Seq1
                        XX = Seq2
                        XX = Seq3
                        X = X
                    End If
                    
                   
                    XOverList(ActiveSeq, SIP).LHolder = 0
                    If LongWindedFlag = 1 And (SEventNumber > 0 Or CircularFlag = 0) Then
                        If XPosDiff(EN) = LenXoverSeq Then
                            ENX = 1
                        Else
                            ENX = XPosDiff(EN) + 1
                        End If
                        If XPosDiff(BE) = 1 Or XPosDiff(BE) = 0 Then
                            BEx = LenXoverSeq
                        Else
                            BEx = XPosDiff(BE) - 1
                        End If
                        
                        'If BE = 1451 And EN = 2460 Then
                        '    X = X
                        'End If
                        If SEventNumber > 0 Then
                            If EWarn = 0 Then Call CheckEndsVB(10, EWarn, LenXoverSeq, 1, CircularFlag, Seq1, Seq2, Seq3, BE, EN, SeqNum(), XPosDiff(), XDiffpos())
                            If BWarn = 0 Then Call CheckEndsVB(10, BWarn, LenXoverSeq, 0, CircularFlag, Seq1, Seq2, Seq3, BE, EN, SeqNum(), XPosDiff(), XDiffpos())
                        End If
                        '1451,2460
                        If BWarn = 1 And EWarn = 1 Then
                            XOverList(ActiveSeq, SIP).SBPFlag = 3
                        ElseIf BWarn = 1 Then
                            XOverList(ActiveSeq, SIP).SBPFlag = 1
                        ElseIf EWarn = 1 Then
                            XOverList(ActiveSeq, SIP).SBPFlag = 2
                        End If
                    End If
                    GoOn = 1
                    XB = BE
                    XE = EN
                    
                    
                    If SEventNumber = 0 And ShowPlotFlag = 2 And (CLine = "" Or CLine = " ") Then
                        StartPlt(8) = 1
                        Call UpdatePlotB(ActiveSeq, ActiveMajorP, ActiveMinorP, SIP)
                    End If
                End If
                
           
    
            End If
            If FindallFlag = 1 Then
                Call SwapRound(WF, WF2, BE, EN, BE2, EN2, nM, nN, nK, nL, odPValue, dPValue2, MultNeg)
                If MCFlag = 0 Then
                    xpvalue = odPValue * MCCorrection
                    
                    If odPValue >= 1 Then
                        dPValue = 1
                    Else
                        dPValue = 1 - (1 - odPValue) ^ MCCorrection
                    End If
                    If dPValue = 0 Then dPValue = xpvalue
                Else
                    xpvalue = odPValue
                End If
            Else
                Exit For
            End If
        Next
    End If

End Sub
Public Function SubPVal(nM, nN, LenXoverSeq, SP, EP, XPosDiff() As Long, XDiffpos() As Long, XoverSeqNumTS() As Long)

Dim X As Long, MinSeen As Long, MaxSeen As Long

If XPosDiff(SP) > 0 Then
    X = XPosDiff(SP) - 1
Else
    X = LenXoverSeq
End If



MaxSeen = -1000
MinSeen = 1000
If XPosDiff(EP) + 1 > LenXoverSeq Then
    epx = LenXoverSeq
Else
    epx = XPosDiff(EP) + 1
End If
If X = epx Then epx = epx - 1
If epx = 0 Then epx = LenXoverSeq

modx = 0
Do While X <> epx
     If X < 0 Then
        X = LenXoverSeq + X
     End If
     If X > LenXoverSeq Then
         modx = XoverSeqNumTS(LenXoverSeq)
         X = X - LenXoverSeq - 1
        
     End If
    If X < 0 Then Exit Function
    If MaxSeen < XoverSeqNumTS(X) + modx Then
        MaxSeen = XoverSeqNumTS(X) + modx: B = X
    End If
    If MinSeen > XoverSeqNumTS(X) + modx Then
        MinSeen = XoverSeqNumTS(X) + modx: E = X
    End If
    If X = epx Then Exit Do
    X = X + 1
Loop
nK = (MaxSeen - MinSeen)

GetTSPVal WF, nM, nN, nK, dPValue

'XX = XTable(26, 47, 36)

SubPVal = dPValue
End Function

Public Sub DeactivateScans()
ConsensusProg = 0
For X = 0 To AddNum - 1
    DoScans(0, X) = 0
Next X
End Sub
Public Sub SSXoverC(FindallFlag, WinNum, SeqMap() As Byte, ZPScoreHolder() As Double, ZSScoreHolder() As Double, CorrectP As Double, oSeq As Long, PermSScores() As Long, PermPScores() As Long, SScoreHolder() As Long, PScoreHolder() As Long, TraceSub() As Long, SeqScore3() As Integer, MeanPScore() As Double, SDPScore() As Double, Seq34Conv() As Byte, VRandConv() As Byte, VRandTemplate() As Byte, HRandTemplate() As Long, TakenPos() As Byte, DG1() As Byte, DG2() As Byte, DoGroupS() As Byte, DoGroupP() As Byte)


Dim WinNumX As Long


WinNumX = WinNum
'61 seconds
Dim NxtX As Long
If SelGrpFlag = 1 Then
    If GrpMaskSeq(Seq1) = 0 And GrpMaskSeq(Seq2) = 0 And GrpMaskSeq(Seq3) = 0 Then
        Exit Sub
    End If
End If

Dim SIP As Long, HVX As Byte, GoOnDraw As Byte, DA As Long, Ma As Long, Mi As Long, maxz As Double, wps As Byte, LSeq As Long, RndNum As Long

Dim DistanceX(2) As Double, ValidX As Double
Dim DoneThis As Byte
Dim SHPos As Long
Dim YPos As Byte, HN1 As Byte, HN2 As Byte, HN3 As Byte, Hi As Integer, Tally() As Long
Dim HV As Integer, LV1 As Integer, LV2 As Integer, LV3 As Integer, HP As Integer, LP1 As Integer, LP2 As Integer
Dim FHPos As Long
Dim tZPScore() As Double, tZSScore() As Double, HRandTemplate2() As Long, LRegion As Long
Dim winscore  As Double, winp As Double, TotP As Double, SP As Long, EP As Long



Rnd (-SSRndSeed)
LSeq = Len(StrainSeq(0))


C = 0

'Set random number seed

'For Seq1 = 0 To NextNo - 2
'    For Seq2 = Seq1 + 1 To NextNo - 1
'        For Seq3 = Seq2 + 1 To NextNo
GoOnDraw = 0
'1.297


Call OrderSeqs(tseq1, tseq2, tseq3, Seq1, Seq2, Seq3, TraceSub())

'1.516
'0.984



If (MaskSeq(Seq1) = 0 And MaskSeq(Seq2) = 0 And MaskSeq(Seq3) = 0 And Seq3 > Seq2 And Seq2 > Seq1) Or (IndividualB > -1 And Seq1 <> Seq2 And Seq2 <> Seq3 And Seq1 <> Seq3 And (Seq1 = IndividualB Or Seq1 = IndividualA) And (Seq2 = IndividualB Or Seq2 = IndividualA) And MaskSeq(Seq3) <= 1) Or (IndividualA > -1 And IndividualB = -1 And Seq1 <> Seq2 And Seq3 > Seq2 And Seq3 <> Seq1 And (Seq1 = IndividualA Or Seq2 = IndividualA) And MaskSeq(Seq3) <= 1 And MaskSeq(Seq2) <= 1) Then
        ReDim Tally(5)
        ReDim PScores(15)
        ReDim tZPScore(1, 15), tZSScore(1, 14)
 '1.687
 '1.484


        WinNumX = CLng(Len(StrainSeq(0)) / SSStep + 1)
        'Dim SeqMap() As Byte, ZPScoreHolder() As Double, ZSScoreHolder() As Double
        
        
        'ReDim SeqMap(Len(StrainSeq(0)))
        'ReDim ZPScoreHolder(WinNumX, 15)
        'ReDim ZSScoreHolder(WinNumX, 14)

        
        Dummy = BlankSSArrays(UBound(ZSScoreHolder, 1), UBound(ZSScoreHolder, 2), UBound(ZPScoreHolder, 1), UBound(ZPScoreHolder, 2), ZSScoreHolder(0, 0), ZPScoreHolder(0, 0))
 
                'Triplet specific stuff

                
                 
                DoneThis = 0
                FHPos = 0
                
                '2.422('8.844/2.375/9.875)?
                '2.172

                'Get Outlyer
                tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                SSOutlyerFlag = pSSOutlyerFlag
                If X = X Then
                    'If UBound(TreeDistance, 1) <> UBound(Distance, 1) Then
                    '    X = X
                    'End If
                    NN1 = UBound(TreeDistance, 1)
                    NN2 = UBound(Distance, 1)
                    'If NN1 <> UBound(TreeDistance, 2) Then
                    '    X = X
                    'End If
                    'If NN2 <> UBound(Distance, 2) Then
                    '    X = X
                    'End If
                    'If PermNextNo > NN1 Or PermNextNo > NN2 Or NextNo > NN1 Or NextNo > NN2 Then
                    '    X = X
                    'End If
                    Dummy = GetSSOL(Len(StrainSeq(0)), SEventNumber, NN1, NN2, PermNextNo, Nextno, SSOutlyerFlag, 0, Seq1, Seq2, Seq3, oSeq, TreeDistance(0, 0), Distance(0, 0), MissingData(0, 0), TraceSub(0))
                Else
                    Call GetSSOutlyer(oSeq, 0, TraceSub())
                End If
                XX = SEventNumber
                If oSeq = -1 Then Exit Sub
                If Seq3 > Nextno Then
                    Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                    SSOutlyerFlag = 0
                    oSeq = Nextno + 1
                End If
                
                '9.156 (12.407?)
                '6.609
                '8.657
                '2.515
  

                
                '@
                
                ValidX = MakeDistX(Seq1, Seq2, Seq3, Len(StrainSeq(0)), DistanceX(0), SeqNum(0, 0))
                    
                '10.781
                '4.062
      
                    
                If FindallFlag = 0 Then
                    If DistanceX(0) < DistanceX(1) And DistanceX(0) < DistanceX(2) Then
                        HV = 0: LV1 = 1: LV2 = 2: LV3 = 3
                        HP = 2: LP1 = 3: LP2 = 4
                    ElseIf DistanceX(1) < DistanceX(0) And DistanceX(1) < DistanceX(2) Then
                        HV = 1: LV1 = 0: LV2 = 2: LV3 = 3
                        HP = 3: LP1 = 2: LP2 = 4
                    Else
                        HV = 2: LV1 = 0: LV2 = 1: LV3 = 3
                        HP = 4: LP1 = 2: LP2 = 3
                    End If
                Else
                    If DistanceX(0) < DistanceX(1) And DistanceX(0) < DistanceX(2) Then
                        LV3 = 0: LV1 = 1: LV2 = 2: HV = 3
                        HP = 2: LP1 = 3: LP2 = 4
                    ElseIf DistanceX(1) < DistanceX(0) And DistanceX(1) < DistanceX(2) Then
                        LV3 = 1: LV1 = 0: LV2 = 2: HV = 3
                        HP = 3: LP1 = 2: LP2 = 4
                    Else
                        LV3 = 2: LV1 = 0: LV2 = 1: HV = 3
                        HP = 4: LP1 = 2: LP2 = 3
                    End If
                End If
                'Calculate 3 sequence scores
                ReDim PScores(Len(StrainSeq(0)))
                
                '10.663
                '4.328
                
                'This could be speeded up with better nesting
                '@'@
                Get3Score SSGapFlag, Len(StrainSeq(0)), Seq1, Seq2, Seq3, SeqNum(0, 0), SeqScore3(0)
                
                '12.047
                '5.844
                
                'Examine the windows (10/18)
                If SSFastFlag = 0 Or ShowPlotFlag = 1 Then
                    GoOn = 1
                Else
                    GoOn = 0
                End If
                For X = 1 To LSeq - SSWinLen Step SSStep
                    
                    
                    If SSFastFlag = 1 And ShowPlotFlag <> 1 Then
                        'Quick check to see if this window is worth examining
                        '@'@
                        NxtX = QuickCheckB(HP, LP1, LP2, SSStep, Len(StrainSeq(0)), X, SSWinLen, Tally(0), SeqScore3(0))
                        If NxtX > LSeq - SSWinLen Then Exit For
                        '@
                        X = NxtX
                        GoOn = 1
                    End If
                    
                    
                     
                     SHPos = X / SSStep
                     '141,161,181,201,261,281,301,321,341,361,381,401,421,441,461,481,501
                     If FHPos = 0 Then FHPos = SHPos
                      
                     If SSOutlyerFlag = 0 Then
                         ReDim PScores(15)
                         
                        GetPScoresRnd X, SSWinLen, Seq1, Seq2, Seq3, oSeq, Len(StrainSeq(0)), PScores(0), Seq34Conv(0, 0), SeqScore3(0), SeqNum(0, 0), HRandTemplate(0)
                     Else
                         
                         
                         '@
                         GetPScores2 X, SSWinLen, Seq1, Seq2, Seq3, oSeq, Len(StrainSeq(0)), PScores(0), Seq34Conv(0, 0), SeqScore3(0), SeqNum(0, 0)
                     End If
                     
                     'Do Perms
                     
                     
                     
                     '@'@'@'@'@
                     Dummy = DoPerms3(Len(StrainSeq(0)), SSWinLen, SSNumPerms, SSNumPerms2, PScores(0), VRandTemplate(0, 0), VRandConv(0, 0), PermPScores(0, 0))
                     
                     'Calculate Z scores
       '
                    '@'@'@'@'@'@
                     Dummy = MakeZValue2(SHPos, WinNumX, 15, SSNumPerms, SSNumPerms2, DG1(0), PermPScores(0, 0), ZPScoreHolder(0, 0))
                     DoSums SSNumPerms, SSNumPerms2, PermSScores(0, 0), PermPScores(0, 0)
                     'Calculate Z Scores for sums
                     '@'@'@'@
                     MakeZValue2 SHPos, WinNumX, 12, SSNumPerms, SSNumPerms2, DG2(0), PermSScores(0, 0), ZSScoreHolder(0, 0)
                     
                    X = X
                Next 'X
                'Exit Sub
                'Exit Sub
                '175.953 : 160.562, 161.125, 161.016, 155.000, 102.093,68.437, 67.531,53.013, 48.738
                '144.313 - MakeZValue2(31.828), 60.641 (7.796), 49.641
                '125.125 - dosums (19.188),52.359 (8.282),41.422
                '87.953 - makexvalue1 (37.172),37.000 (15.359),37.203
                '47.672 - doperms3 (40.281),19.329 (17.671) 19.266
                '35.188 - getpscores2 (12.484),14.109 (5.22)
                '34.063- quickcheck (1.125)
                'Background (22.012):9.000
                
         
                'GoOn = 0
                If GoOn = 1 Then
                    GoOn = 0
                    
                    
                    'Find potentially recombinant regions from the S and P plots
                    
                    
                    
                    If FindallFlag = 0 Then
                        Dummy = FindMaxZ(HV, LV1, LV2, LV3, WinNumX, FHPos, SHPos, SeqMap(0), DoGroupS(0, 0), DoGroupP(0, 0), maxz, winp, winscore, wps, ZSScoreHolder(0, 0), ZPScoreHolder(0, 0))
                    Else
                        Dummy = FindMaxZ(HV, LV1, LV2, LV3, WinNumX, FHPos, SHPos, SeqMap(0), DoGroupS(0, 0), DoGroupP(0, 0), maxz, winp, winscore, wps, ZSScoreHolder(0, 0), ZPScoreHolder(0, 0))
                    End If
                    
                                    
                    'Find the regions (ie not just the interesting windows)
                    'This is done by expanding windows until P-values stop increasing.
                    For X = FHPos To SHPos
                        
                        If SeqMap(X) <> HV Then
                            If SeqMap(X) <> 3 Then
                                HVX = SeqMap(X)
                                
                                EP = 1
                                SP = X
                                Do
                                    '@
                                    If SeqMap(X + EP) <> SeqMap(SP) Then
                                        If FindallFlag = 0 Then
                                            Exit Do
                                        Else
                                            If SeqMap(X + EP) <> HV Then
                                                Exit Do
                                            End If
                                        End If
                                    End If
                                    EP = EP + 1
                                    If X + EP > SHPos Then
                                        If CircularFlag = 0 Then
                                            EP = SHPos - X + 1
                                            Exit Do
                                        Else
                                            If SHPos = WinNumX Then
                                                EP = 1
                                                X = FHPos
                                            Else
                                                EP = SHPos - X + 1
                                                Exit Do
                                            End If
                                        End If
                                    End If
                                Loop
                                'XX = (EP * 20)
                                EP = X + EP - 1
                                'Erase trace of region so it isn't redone
                                If EP >= SP Then
                                    For z = SP To EP
                                        SeqMap(z) = HV
                                    Next z
                                Else
                                    For z = 0 To EP
                                        SeqMap(z) = HV
                                    Next z
                                    For z = EP To WinNumX
                                        SeqMap(z) = HV
                                    Next z
                                End If
                                
                                'Shrink region
                                
                                
                                Dim TBegin As Long, TEnd As Long
                                TBegin = 0
                                TEnd = 0
                                'HVX = SeqMap(EP)
                                'B = TBegin
                                
                                '4:39
                                '835 : 4:03 : 4:01
                                If X = X Then
                                     If HVX = 3 Then
                                        X = X
                                    End If
                                    Dummy = ShrinkRegionC(Len(StrainSeq(0)), SSGapFlag, SSStep, SSWinLen, HVX, Seq1, Seq2, Seq3, EP, SP, TEnd, TBegin, SeqNum(0, 0))
                                Else
                                
                                    Call ShrinkRegion(HVX, TBegin, TEnd, SP, EP, SSStep, SSWinLen, SSGapFlag, SeqNum())
                                End If
                                If TBegin <> 0 And TEnd <> 0 And TBegin <> TEnd Then
                                    'Work out significance
                                    
                                    If TBegin < TEnd Then
                                        LRegion = TEnd - TBegin + 1
                                    Else
                                        LRegion = TEnd + (Len(StrainSeq(0)) - TBegin + 1)
                                    End If
                                    
                                    If SSOutlyerFlag = 0 Then
                                        ReDim PScores(15), HRandTemplate2(LRegion), TakenPos(LRegion)
                                        For z = 1 To LRegion
                                            RndNum = Int((LRegion * Rnd) + 1)
                                            If TakenPos(RndNum) = 0 Then
                                                HRandTemplate2(z) = RndNum
                                                TakenPos(RndNum) = 1
                                            Else 'find next available position to the right
                                                Y = RndNum
                                                Do While TakenPos(Y) = 1
                                                    Y = Y + 1
                                                    If Y > LRegion Then Y = 1
                                                Loop
                                                HRandTemplate2(z) = Y
                                                TakenPos(Y) = 1
                                            End If
                                        Next 'X
                                        
                                        GetPScoresRnd TBegin, LRegion, Seq1, Seq2, Seq3, oSeq, Len(StrainSeq(0)), PScores(0), Seq34Conv(0, 0), SeqScore3(0), SeqNum(0, 0), HRandTemplate2(0)
                                                    
                                    Else
                                        
                                        
                                        GetPScores2 TBegin, LRegion, Seq1, Seq2, Seq3, oSeq, Len(StrainSeq(0)), PScores(0), Seq34Conv(0, 0), SeqScore3(0), SeqNum(0, 0)
                                    End If
                                                
                                   XX = SEventNumber
                                    
                                    '@'@'@'@'@'@'@'@'@'@'@
                                    Dummy = DoPerms3(Len(StrainSeq(0)), LRegion, SSNumPerms, SSNumPerms, PScores(0), VRandTemplate(0, 0), VRandConv(0, 0), PermPScores(0, 0))
                                    '@'@
                                    MakeZValue2 0, 1, 15, SSNumPerms, SSNumPerms, DG1(0), PermPScores(0, 0), tZPScore(0, 0)
                                    
                                    '6:25 - 200 win 200,20 perms
                                    '11:12 - 200 win 1000,100 perms
                                    
                                    
                                    DoSums SSNumPerms, SSNumPerms, PermSScores(0, 0), PermPScores(0, 0)
                                        
                                    'SScoreHolder(SHPos, 1) = PermSScores(0, 1)
                                    'SScoreHolder(SHPos, 2) = PermSScores(0, 2)
                                    'SScoreHolder(SHPos, 3) = PermSScores(0, 3)
                                    'SScoreHolder(SHPos, 4) = PermSScores(0, 4)
                                    'SScoreHolder(SHPos, 5) = PermSScores(0, 5)
                                    'SScoreHolder(SHPos, 7) = PermSScores(0, 7)
                                            
                                    'Calculate Z Scores
                                    '@'@'@
                                    MakeZValue2 0, 1, 12, SSNumPerms, SSNumPerms, DG2(0), PermSScores(0, 0), tZSScore(0, 0)
                                                
                                    maxz = 0: winp = 0: wps = 0
                                                
                                    For z = 0 To 1
                                        If Abs(tZPScore(0, DoGroupP(z, HV))) < Abs(tZPScore(0, DoGroupP(z, LV1))) Or Abs(tZPScore(0, DoGroupP(z, HV))) < Abs(tZPScore(0, DoGroupP(z, LV2))) Then
                                                        
                                            For Y = 0 To 2
                                                If Y <> HV Then
                                                    If Abs(tZPScore(0, DoGroupP(z, Y))) > maxz Then
                                                        maxz = Abs(tZPScore(0, DoGroupP(z, Y)))
                                                        winscore = DoGroupP(z, Y)
                                                        wps = 1
                                                    End If
                                                End If
                                            Next Y
                                        End If
                                    Next z
                                    
                                    For z = 0 To 1
                                        If Abs(tZSScore(0, DoGroupS(z, HV))) < Abs(tZSScore(0, DoGroupS(z, LV1))) Or Abs(tZSScore(0, DoGroupS(z, HV))) < Abs(tZSScore(0, DoGroupS(z, LV2))) Then
                                            For Y = 0 To 2
                                                If Y <> HV Then
                                                    If Abs(tZSScore(0, DoGroupS(z, Y))) > maxz Then
                                                        maxz = Abs(tZSScore(0, DoGroupS(z, Y)))
                                                        winscore = DoGroupS(z, Y)
                                                        wps = 2
                                                    End If
                                                End If
                                            Next Y
                                        End If
                                    Next z
                                            
                                    If maxz > CriticalZ Or (ShortOutFlag = 3 And maxz <> 0) Then
                                        
                                        If TBegin < TEnd Then
                                            rlen = TEnd - TBegin
                                        Else
                                            rlen = TEnd + (Len(StrainSeq(0)) - TBegin)
                                        End If
                                        winp = NormalZ(maxz)
                                        winp = winp * (Len(StrainSeq(0)) / rlen)
                                        If MCFlag = 0 Then
                                            winp = winp * MCCorrection * (Len(StrainSeq(0)) / SSWinLen)
                                        End If
                                        If ShortOutFlag = 3 Then
                                            If winp <= mtP(5) Then
                                                mtP(5) = winp
                                            End If
                                        End If
                                        XX = SEventNumber
                                        If winp < LowestProb Then
                                            
                                            GoOnDraw = 1
                                            'Keep track of signal numbers
                                            oRecombNo(100) = oRecombNo(100) + 1
                                            oRecombNo(5) = oRecombNo(5) + 1
                                            If ShortOutFlag = 1 Then
                                                ShortOutput(5) = 1
                                                AbortFlag = 1
                                                Exit Sub
                                            End If
                                            
                                            ReDim Tally(5)
                                            If TBegin < TEnd Then
                                                QuickCheck TBegin, TEnd - TBegin, Tally(0), SeqScore3(0)
                                            Else
                                                QuickCheck 1, TEnd, Tally(0), SeqScore3(0)
                                                QuickCheck TBegin, Len(StrainSeq(0)), Tally(0), SeqScore3(0)
                                            End If
                                                        
                                            If HP = 2 Then
                                                If Tally(LP1) > Tally(LP2) Then
                                                    ActiveSeq = Seq1: ActiveMajorP = Seq2: ActiveMinorP = Seq3
                                                Else
                                                    ActiveSeq = Seq2: ActiveMajorP = Seq1: ActiveMinorP = Seq3
                                                End If
                                            ElseIf HP = 3 Then
                                                If Tally(LP1) > Tally(LP2) Then
                                                    ActiveSeq = Seq1: ActiveMajorP = Seq3: ActiveMinorP = Seq2
                                                Else
                                                    ActiveSeq = Seq3: ActiveMajorP = Seq1: ActiveMinorP = Seq2
                                                End If
                                            ElseIf HP = 4 Then
                                                If Tally(LP1) > Tally(LP2) Then
                                                    ActiveSeq = Seq2: ActiveMajorP = Seq3: ActiveMinorP = Seq1
                                                Else
                                                    ActiveSeq = Seq3: ActiveMajorP = Seq2: ActiveMinorP = Seq1
                                                End If
                                            End If
                                            SIP = 0
                                            If APermFlag = 0 Then
                                                Call UpdateXOList3(ActiveSeq, CurrentXOver(), XOverList(), 5, winp, SIP)
                                            Else
                                                SIP = 1
                                            End If
                                            'Call UpdateXOList(ActiveSeq, CurrentXOver(), XoverList())
                                            If SIP > 0 Then
                                                XOverList(ActiveSeq, SIP).LHolder = winscore + 10 * (wps - 1)
                                                XOverList(ActiveSeq, SIP).MajorP = ActiveMajorP
                                                XOverList(ActiveSeq, SIP).MinorP = ActiveMinorP
                                                XOverList(ActiveSeq, SIP).DHolder = oSeq
                                                XOverList(ActiveSeq, SIP).Daughter = ActiveSeq
                                                LenXoverSeq = BSSubSeq(Len(StrainSeq(0)), Seq1, Seq2, Seq3, SeqNum(0, 0), XPosDiff(0), XDiffpos(0), Scores(0, 0))
                                                
                                                
                                                
                                                EWarn = 0: BWarn = 0
                                                
                                                Call CentreBP(TBegin, TEnd, XPosDiff(), XDiffpos(), BWarn, EWarn, 10, LenXoverSeq)
                                                
                                                XOverList(ActiveSeq, SIP).Beginning = TBegin
                                                XOverList(ActiveSeq, SIP).Ending = TEnd
                                                
                                                Call FixEnds(BWarn, EWarn, MissingData(), XOverList(), ActiveSeq, SIP)
                                                
                                                'XOverList(ActiveSeq, SIP).Beginning = TBegin
                                                'XOverList(ActiveSeq, SIP).Ending = TEnd
                                                
                                                XOverList(ActiveSeq, SIP).Probability = winp
                                                XOverList(ActiveSeq, SIP).ProgramFlag = 5
                                                
                                                If TBegin = 478 And TEnd = 994 Then
                                                    X = X
                                                End If
                                                If MCFlag = 2 Then
                                                    ProbabilityXOver = winp
                                                    If -Log10(ProbabilityXOver) * 2 > 0 And -Log10(ProbabilityXOver) * 2 < 100 Then
                                                        PValCat(CurrentCorrect, CInt(-Log10(ProbabilityXOver) * 2)) = PValCat(CurrentCorrect, CInt(-Log10(ProbabilityXOver) * 2)) + 1
                                                    ElseIf CInt(-Log10(ProbabilityXOver) * 2) >= 100 Then
                                                        PValCat(CurrentCorrect, 100) = PValCat(CurrentCorrect, 100) + 1
                                                    End If
                                                End If
                                                If SEventNumber = 0 And ShowPlotFlag = 2 And (CLine = "" Or CLine = " ") Then
                                                    StartPlt(5) = 1
                                                    Call UpdatePlotB(ActiveSeq, ActiveMajorP, ActiveMinorP, SIP)
                                                    
                                                End If
                                            End If
                                            
                                            
                                            
                                            'Call FindDaughter(ActiveSeq, ActiveMinorP, ActiveMajorP, 0, 0, 5, SIP)
                                                        
                                        End If
                                    End If
                                End If
                            End If
                        End If
                    Next X
                End If
                If ShowPlotFlag = 1 And GoOnDraw = 1 Then
                    
                    GoOnDraw = 0
                    
                                        
                'findmax and min
                    
                                            
                    Call DrawSSPlots(DG1(), ZPScoreHolder(), ZSScoreHolder())
                                          
                End If
Skip:
                C = C + 1
            End If
            
'            ET = GetTickCount()
'
'            If ET - LT > 500 Or AbortFlag = 1 Then
'                SSE = GetTickCount()
'                If C < MCCorrection + 1 Then
'                    Form1.ProgressBar1.Value = (C + 1) / (MCCorrection + 1) * 100
'                End If
'                Form1.Label69(0).Caption = DoTimeII(SSE - SSS)
'                Form1.Label57(0).Caption = DoTimeII(ET - ST)
'                Form1.Label69(1).Caption = CStr(TotalSSRecombinants)
'                Form1.Label57(1).Caption = CStr( oRecombNo(100) )
'                If oRecombNo(100) > oRec And ShowPlotFlag = 2 And (CLine = "" Or CLine = " ") Then
'                    StartPlt(5) = 1
 '                   oRec =  oRecombNo(100)
'                    Call UpdatePlotC
'                End If
'                LT = ET
'                Form1.SSPanel13.Caption = "Approximately " & DoTime((SSE - SSS) * (100 / Form1.ProgressBar1.Value) - (SSE - SSS)) & " remaining"
'                Form1.SSPanel1.Caption = Str(C) & " of" & Str(MCCorrection) & " triplets examined"
'                DoEvents
'
'                If AbortFlag = 1 Then
'                    Form1.SSPanel1.Caption = ""
'                    Form1.ProgressBar1.Value = 0
'                    Exit Sub
'                End If
'
'            End If
'
ESub:
            
            Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
'        Next 'Seq3
'    Next 'Seq2
'Next 'Seq1


End Sub
Private Function GetHTMLFromURL(sUrl As String) As String
'From Peter G. Aitken
Dim s As String
Dim hOpen As Long
Dim hOpenUrl As Long
Dim bDoLoop As Boolean
Dim bRet As Boolean
Dim sReadBuffer As String * 2048
Dim lNumberOfBytesRead As Long

hOpen = InternetOpen(scUserAgent, INTERNET_OPEN_TYPE_PRECONFIG, vbNullString, vbNullString, 0)
hOpenUrl = InternetOpenUrl(hOpen, sUrl, vbNullString, 0, INTERNET_FLAG_RELOAD, 0)

bDoLoop = True
While bDoLoop
    sReadBuffer = vbNullString
    bRet = InternetReadFile(hOpenUrl, sReadBuffer, Len(sReadBuffer), lNumberOfBytesRead)
    s = s & left$(sReadBuffer, lNumberOfBytesRead)
    If Not CBool(lNumberOfBytesRead) Then bDoLoop = False
Wend

If hOpenUrl <> 0 Then InternetCloseHandle (hOpenUrl)
If hOpen <> 0 Then InternetCloseHandle (hOpen)

GetHTMLFromURL = s

End Function


Public Sub NucleotideEg()


Dim CNum As Long, WebEnv As String, QueryKey As String, RetMax As String, lString As String, Utils As String, db As String, Query As String, Report As String, ESearch As String, Efetch As String

Utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils"
db = "nucleotide"
Query = "begomovirus"
'docsum, brief, gb, citation, medline, asn.1, mlasn1, uilist, sgml, gen, fasta
Report = "fasta"
RetMax = "20"
ESearch = Utils + "/esearch.fcgi?&email=darren@science.uct.ac.za&tool=darrenmartinrdp&db=" + db + "&usehistory=y&term=" + Query
lString = GetHTMLFromURL(ESearch) 'Form1.Inet1.OpenURL(ESearch)
If Len(lString) > 0 Then
    Pos = InStr(1, lString, "<Count>", vbBinaryCompare)
    Pos2 = InStr(1, lString, "</Count>", vbBinaryCompare)
    CNum = Mid$(lString, Pos + 7, (Pos2 - (Pos + 7)))
    Pos = InStr(1, lString, "<QueryKey>", vbBinaryCompare)
    Pos2 = InStr(1, lString, "</QueryKey>", vbBinaryCompare)
    QueryKey = Mid$(lString, Pos + 10, (Pos2 - (Pos + 10)))
    Pos = InStr(1, lString, "<WebEnv>", vbBinaryCompare)
    Pos2 = InStr(1, lString, "</WebEnv>", vbBinaryCompare)
    WebEnv = Mid$(lString, Pos + 8, Pos2 - (Pos + 8))
    'lString = ""
    ChDir App.Path
    On Error Resume Next
    Kill Query + ".fas"
    On Error GoTo 0
    Open Query + ".fas" For Output As #1
    
    'For X = 0 To CNum - 1
    '    Efetch = Utils + "/efetch.fcgi?&email=darren@science.uct.ac.za&tool=darrenmartinrdp&rettype=" + Report + "&retmode=text&retstart=" + Trim(Str(X)) + "&retmax=1&db=" + db + "&query_key=" + QueryKey + "&WebEnv=" + WebEnv
    '    lString = Form1.Inet1.OpenURL(Efetch)
    '
    '    Put #1, , lString
    'Next X
    'get in 1 chunk
    'For X = 0 To CNum - 1
    lString = ""
    
    Dim Bitesize() As Byte, LenToss As Long, RS As Long, Target As String, T2 As String, lString2 As String, LoadChunk As Long
    
    For X = 0 To CNum - 1 Step RetMax
        'get the sequence title
        ST = 1: EN = 2
        Efetch = Utils + "/efetch.fcgi?&email=darren@science.uct.ac.za&tool=darrenmartinrdp&retmax=" + RetMax + "&retmode=text&retstart=" + Trim(Str(X)) + "&rettype=" + Report + "&db=" + db + "&query_key=" + QueryKey + "&WebEnv=" + WebEnv
        lString = ""
        z = 0
        Do While lString = "" Or left(lString, 5) = "Error"
            lString = GetHTMLFromURL(Efetch) 'Form1.Inet1.OpenURL(Efetch)
            z = z + 1
            If z = 20 Then
                Exit Sub
                XX = Len(lString)
            End If
        Loop
        Print #1, lString
        
        
        If X + 3 < CNum Then
            Form1.ProgressBar1.Value = (X + 3) / CNum * 100
            Form1.SSPanel1.Caption = (X + 3) & " of " & CNum & " sequences loaded"
        Else
            Form1.ProgressBar1.Value = 100
            Form1.SSPanel1.Caption = CNum & " of " & CNum & " sequences loaded"
        End If
    Next X
    XX = Len(lString)
     
     
    'Next X
    Close #1
End If
Form1.SSPanel1.Caption = ""
Form1.ProgressBar1 = 0
'Open "testout" For Output As #1
'Print #1, Efetch
'Print #1, lString
'Close #1
X = X
End Sub

Public Sub SaveAlign(SIndex As Integer)
    'Subroutine used to save multiple sequence alignments
    'I must remember to incorporate the "save recombinant region" thing here
    Dim SeqSave() As Integer
    Dim SeqMask() As Byte, DontInclude As Byte
    Dim Addj As Integer, Addon As Integer
    Dim NumSeqs As Long, LastPos As Long, X As Long, Y As Long, z As Long
    Dim TString As String, AName As String
    Dim Namestring() As String, TempSeq() As String, TempName() As String

    ReDim TempSeq(Nextno)
    ReDim TempName(Nextno)

    With Form1.CommonDialog1
        .FileName = ""
        .DefaultExt = "fas"
        .Filter = "All Files (*.*)|*.*|RDP Project Files (*.rdp)|*.rdp|Clustal Multiple Alignment Format (*.aln)|*.aln|DNAMan Multiple Alignment Format (*.msd)|*.msd|FASTA Multiple Alignment Format (*.fas)|*.fas|GCG Multiple Alignment Format (*.gcg)|*.gcg|GDE Multiple Alignment Format (*.gde)|*.gde|NEXUS Format (*.nex)|*.nex|Mega Multiple Alignment Format (*.meg)|*.meg|PHYLIP Format (*.phy)|*.phy|Sequence Files (*.seq)|*.seq"
        .Action = 2 'Specify that the "open file" action is required.
        AName = .FileName  'Stores selected file name in the
        ANameII = .FileTitle
    End With

    If AName = "" Then Exit Sub
    NumSeqs = 0

    If SIndex = 0 Then
        NumSeqs = Nextno
        'zz = 0
        For X = 0 To Nextno
            
            TempName(X) = StraiName(X)
            TempSeq(X) = StrainSeq(X)
           
        Next 'X
    ElseIf SIndex = 5 Then ' save only disabled
        For X = 0 To Nextno

            If MaskSeq(X) > 0 Then
                TempSeq(NumSeqs) = StrainSeq(X)
                TempName(NumSeqs) = StraiName(X)
                NumSeqs = NumSeqs + 1
            End If

        Next 'X

        NumSeqs = NumSeqs - 1
    ElseIf SIndex = 6 Then ' save only enabled
        For X = 0 To Nextno

            If MaskSeq(X) = 0 Then
                TempSeq(NumSeqs) = StrainSeq(X)
                TempName(NumSeqs) = StraiName(X)
                NumSeqs = NumSeqs + 1
            End If

        Next 'X

        NumSeqs = NumSeqs - 1
    ElseIf SIndex = 1 Then

        For X = 0 To Nextno

            If CurrentXOver(X) = 0 Then
                TempSeq(NumSeqs) = StrainSeq(X)
                TempName(NumSeqs) = StraiName(X)
                NumSeqs = NumSeqs + 1
            Else
                DontInclude = 0

                For Y = 1 To CurrentXOver(X)
                    If XOverList(X, Y).Accept = 1 Then
                    
                        DontInclude = 1
                        Exit For
                    End If

                Next 'Y

                If DontInclude = 0 Then
                    TempSeq(NumSeqs) = StrainSeq(X)
                    TempName(NumSeqs) = StraiName(X)
                    NumSeqs = NumSeqs + 1
                End If

            End If

        Next 'X

        NumSeqs = NumSeqs - 1
    ElseIf SIndex = 2 Then
        NumSeqs = Nextno
        ReDim SeqSave(Len(StrainSeq(0)), Nextno)
        For X = 0 To Nextno
            For Y = 1 To Len(StrainSeq(0))
                SeqSave(Y, X) = PermSeqNum(Y, X)
            Next Y
        Next X
        For X = 0 To Nextno
            For Y = 1 To CurrentXOver(X)
                If XOverList(X, Y).Accept = 1 Then
                    BE = XOverList(X, Y).Beginning
                    EN = XOverList(X, Y).Ending
                    
                    For A = 0 To Nextno
                        If BE < EN Then
                            For z = BE To EN
                                SeqSave(z, A) = 46
                            Next z
                        Else
                            For z = BE To Len(StrainSeq(0))
                                SeqSave(z, A) = 46
                            Next z
                            For z = 1 To EN
                                SeqSave(z, A) = 46
                            Next z
                        End If
                    Next A
                End If
            Next Y
        
        Next X
        For X = 0 To Nextno
            TempSeq(X) = StrainSeq(X)
            TempName(X) = StraiName(X)
            For Y = 1 To Len(StrainSeq(0))
                If SeqSave(Y, X) <> PermSeqNum(Y, X) Then
                    Mid$(TempSeq(X), Y, 1) = "-"
                End If
            Next Y
        Next X

    ElseIf SIndex = 3 Then
        NumSeqs = Nextno
        
        
        ReDim SeqSave(Len(StrainSeq(0)), Nextno)
        For X = 0 To Nextno
            For Y = 1 To Len(StrainSeq(0))
                SeqSave(Y, X) = PermSeqNum(Y, X)
            Next Y
        Next X
        For X = 0 To Nextno
            For Y = 1 To CurrentXOver(X)
                If XOverList(X, Y).Accept = 1 Then
                    BE = XOverList(X, Y).Beginning
                    EN = XOverList(X, Y).Ending
                    If BE < EN Then
                        For z = BE To EN
                            SeqSave(z, X) = 46
                        Next z
                    Else
                        For z = BE To Len(StrainSeq(0))
                            SeqSave(z, X) = 46
                        Next z
                        For z = 1 To EN
                            SeqSave(z, X) = 46
                        Next z
                    End If
                End If
            Next Y
        
        Next X
        For X = 0 To Nextno
            TempSeq(X) = StrainSeq(X)
            TempName(X) = StraiName(X)
            For Y = 1 To Len(StrainSeq(0))
                If SeqSave(Y, X) <> PermSeqNum(Y, X) Then
                    Mid$(TempSeq(X), Y, 1) = "-"
                End If
            Next Y
        Next X
        
    ElseIf SIndex = 7 Then
        NumSeqs = Nextno
        
        
        ReDim SeqSave(Len(StrainSeq(0)), Nextno)
        For X = 0 To Nextno
            For Y = 1 To Len(StrainSeq(0))
                SeqSave(Y, X) = PermSeqNum(Y, X)
            Next Y
        Next X
        For X = 0 To Nextno
            For Y = 1 To CurrentXOver(X)
                If XOverList(X, Y).Accept = 1 Then
                    Nextno = Nextno + 1
                    
                    ReDim Preserve StraiName(Nextno), SeqSave(Len(StrainSeq(0)), Nextno)
                    StraiName(Nextno) = StraiName(X) + Str(Nextno)
                    BE = XOverList(X, Y).Beginning
                    EN = XOverList(X, Y).Ending
                    
                    If BE < EN Then
                        For z = BE To EN
                            SeqSave(z, Nextno) = SeqSave(z, X)
                            SeqSave(z, X) = 46
                        Next z
                    Else
                        For z = BE To Len(StrainSeq(0))
                            SeqSave(z, Nextno) = SeqSave(z, X)
                            SeqSave(z, X) = 46
                        Next z
                        For z = 1 To EN
                            SeqSave(z, Nextno) = SeqSave(z, X)
                            SeqSave(z, X) = 46
                        Next z
                    End If
                End If
            Next Y
        
        Next X
        ReDim TempSeq(Nextno)
        ReDim TempName(Nextno)
        For X = 0 To Nextno
            TempSeq(X) = ""
            TempName(X) = StraiName(X)
            For Y = 1 To Len(StrainSeq(0))
                If SeqSave(Y, X) = 0 Then SeqSave(Y, X) = 46
                TempSeq(X) = TempSeq(X) + Chr(SeqSave(Y, X) - 1)
            Next Y
        Next X
        X = 0
        Do While X < Nextno
            If TempSeq(X) = String(Len(StrainSeq(0)), "-") Then
                If X < Nextno Then
                    TempSeq(X) = TempSeq(Nextno)
                    StraiName(X) = StraiName(Nextno)
                End If
                Nextno = Nextno - 1
                X = X - 1
            End If
            X = X + 1
        Loop
        NumSeqs = Nextno
    End If

    
    Close #1
    
    Open AName For Output As #1
    If right$(AName, 4) = ".msd" Or right$(AName, 4) = ".MSD" Then
        'Save DNAMAN file
        Print #1, "FILE: Multiple_Sequence_Alignment"
        Print #1, "PROJECT:"
        Print #1, "NUMBER:" + CStr(NumSeqs + 1)
        Print #1, "MAXLENGTH:" + CStr(Len(TempSeq(0)))
        Addj = 65

        For X = 0 To NumSeqs

            If Len(TempName(X)) > 9 Then
                TempName(X) = Mid$(TempName(X), 1, 8) + Chr$(Addj)
                Addj = Addj + 1
            End If

        Next 'X

        ReDim Namestring(Int((NumSeqs + 1) / 6) + 1)
        Namestring(0) = "NAMES:"

        For X = 1 To Int((NumSeqs + 1) / 6)
            Namestring(X) = "      "
        Next 'X

        For X = 0 To Int((NumSeqs + 1) / 6)
            Y = 0

            Do Until Y = 6 Or (X * 6 + Y) > Nextno
                Namestring(X) = Namestring(X) + " " + TempName(X * 6 + Y)
                Y = Y + 1
            Loop

        Next 'X

        For X = 0 To Int((Nextno + 1) / 6)
            Print #1, Namestring(X)
        Next 'X

        Print #1, ""
        Print #1, "ORIGIN"

        For X = 0 To Nextno
            Addon = 10 - Len(TempName(X))
            TempName(X) = TempName(X) + String$(Addon, 32)
        Next 'X

        For X = 1 To Int(Len(TempSeq(0)) / 60)

            For Y = 0 To NumSeqs
                Print #1, TempName(Y) + Mid$(TempSeq(Y), X * 60 - 59, 60)
            Next 'Y

            Print #1, ""
        Next 'X

        For Y = 0 To NumSeqs
            Print #1, TempName(Y) + Mid$(TempSeq(Y), X * 60 - 59, Len(TempSeq(0)) - (X * 60 - 60))
        Next 'Y

    ElseIf right$(AName, 4) = ".aln" Or right$(AName, 4) = ".ALN" Then
        Addj = 0

        For X = 0 To NumSeqs

            If Len(TempName(X)) > 15 Then
                TempName(X) = Mid$(TempName(X), 1, 14) + Chr$(Addj)
                Addj = Addj + 1
            End If

            Addon = 16 - Len(TempName(X))
            TempName(X) = TempName(X) + String$(Addon, 32)
        Next 'X

        Print #1, "CLUSTAL multiple sequence alignment"
        Print #1, ""
        Print #1, ""

        For Y = 1 To Len(TempSeq(0)) Step 60

            For X = 0 To NumSeqs
                Print #1, TempName(X) + Mid$(TempSeq(X), Y, 60)
            Next 'X

            Print #1, ""
            Print #1, ""
        Next 'Y

    ElseIf right$(AName, 4) = ".phy" Or right$(AName, 4) = ".PHY" Then
        'Save alignment in phylip Format
        Addj = 65
        maxns = 0
        For X = 0 To NumSeqs
            If maxns < Len(TempName(X)) Then maxns = Len(TempName(X))
        Next X
        maxns = maxns + 1
        maxns = 10
        If X = X Then
            Print #1, Trim$(CStr(NumSeqs + 1)) & " " & Trim$(CStr(Len(TempSeq(0)))) & " 1"
            For X = 0 To NumSeqs
                If Len(TempName(X)) < 10 Then
                    Print #1, TempName(X) + String(maxns - Len(TempName(X)), " ") + TempSeq(X)
                Else
                    Print #1, left(TempName(X), 10) + TempSeq(X)
                End If
            Next 'X
            Close #1
            Exit Sub
        Else
            For X = 0 To NumSeqs
    
                If Len(TempName(X)) > 10 Then
                    TempName(X) = Mid$(TempName(X), 1, 8) + Chr$(Addj)
                    Addj = Addj + 1
                End If
    
                Addon = 10 - Len(TempName(X))
                TempName(X) = TempName(X) + String$(Addon, 32)
            Next 'X
        End If
        Print #1, "    " & Trim$(CStr(NumSeqs + 1)) & "   " & Trim$(CStr(Len(TempSeq(0))))

        For Y = 1 To Len(TempSeq(0)) Step 50

            For X = 0 To NumSeqs
                TString = ""

                For z = 0 To 49 Step 10
                    TString = TString + " " + Mid$(TempSeq(X), Y + z, 10)
                Next 'Z

                Print #1, TempName(X) + TString

                If Y = 1 Then
                    TempName(X) = String$(10, " ")
                End If

            Next 'X

            Print #1, ""
        Next 'Y

    ElseIf right$(AName, 4) = ".gde" Or right$(AName, 4) = ".GDE" Then
        'Save alignment in GDE Format

        For X = 0 To NumSeqs
            Print #1, "#" & TempName(X)

            For Y = 1 To Len(TempSeq(0)) + 60 Step 60
                Print #1, Mid$(TempSeq(X), Y, 60)
            Next 'Y

        Next 'X

    Else
        'Save alignment in FASTA format

        For X = 0 To NumSeqs
            Print #1, ">" & TempName(X)

            For Y = 1 To Len(TempSeq(0)) + 70 Step 70
                Print #1, Mid$(TempSeq(X), Y, 70)
            Next 'Y

        Next 'X

    End If
    If SIndex = 7 Then Call UnModSeqNum(0)

    Close #1
End Sub



Public Sub FindSimilar(DaughterSeq, RecNumber, AcceptFlag, StartAcceptX)

Dim LenXOver2 As Long, LenXOver As Long, HitsPerseq() As Long, RSP As Long, REP As Long, RSP2 As Long, REP2 As Long, RS As Long, RE As Long, TP As Long, RS2 As Long, RE2 As Long, TP2 As Long, LSeq As Long
Dim GroupEvents() As Long, Eventholder As Long, EventPos() As Long, EventTraceA() As Long, EventTraceB() As Long
Dim EventBegin() As Long, EventEnd() As Long, TreePoses() As Long, PosNumber() As Long
Dim CNum() As Long, CPVal() As Double

LSeq = Len(StrainSeq(0))
ReDim GEvents(CurrentXOver(DaughterSeq))
ENumbs = 0
X = DaughterSeq
Y = RecNumber
StartAccept = StartAcceptX
            
            
            RSP = XOverList(X, Y).Beginning
            REP = XOverList(X, Y).Ending
'XXXZZZ              TP = XOverlist(X, Y).TreePos(0)
            
                If RSP > REP Then REP = LSeq + REP
                LenXOver = REP - RSP
                
                RE = REP - LenXOver / 4
                RS = RSP + LenXOver / 4
                
                For z = 1 To CurrentXOver(X)
                    
                    If XOverList(X, z).Accept = StartAccept Then
                        'Exit Sub
                        If X = 123456 Then
                            RSP2 = XOverList(X, z).Beginning
                            REP2 = XOverList(X, z).Ending
                            
'XXXZZZ                                  TP2 = XOverlist(X, Z).TreePos(0)
                                
                                If RSP2 < REP2 Then
                                    If REP > LSeq Then
                                        If REP2 < RSP Then
                                            RSP2 = RSP2 + LSeq
                                            REP2 = REP2 + LSeq
                                        End If
                                    End If
                                Else
                                    If REP > LSeq Then
                                        REP2 = REP2 + LSeq
                                    Else
                                        RSP2 = RSP2 - LSeq
                                    End If
                                End If
                                
                                LenXOver2 = REP2 - RSP2
                                OL = 0
                                ol2 = 0
                                If REP > REP2 Then
                                    OL = REP - REP2
                                ElseIf REP < REP2 Then
                                    ol2 = REP2 - REP
                                End If
                                If RSP > RSP2 Then
                                    ol2 = ol2 + RSP - RSP2
                                ElseIf RSP < RSP2 Then
                                    OL = OL + RSP2 - RSP
                                End If
                                RE2 = REP2 - LenXOver2 / 4
                                RS2 = RSP2 + LenXOver2 / 4
                                'Exit Sub
                                
                                If LenXOver > 0 And LenXOver2 > 0 Then
                                    If RS2 < RE And RE2 > RS And TP = TP2 And OL / LenXOver < 0.5 And ol2 / LenXOver2 < 0.5 Then
                                        ENumbs = ENumbs + 1
                                        GEvents(ENumbs) = z
                                        XOverList(X, z).Accept = AcceptFlag
                                        
                                    End If
                                End If
                            Else
                                If XOverList(X, z).Eventnumber = XOverList(X, Y).Eventnumber Then
                                    ENumbs = ENumbs + 1
                                    GEvents(ENumbs) = z
                                    XOverList(X, z).Accept = AcceptFlag
                                End If
                            End If
                    End If
                Next z
End Sub
Public Sub ReCheck()

    Dim HitsPerseq() As Long
    Dim NumHits() As Long
    Dim LSeq As Long, BPos As Long, EPos As Long, Ma As Long, Mi As Long, SubM As Long
    Dim X As Long, EY As Long, Y As Long, z As Long
    Dim NumToAdd As Long
    NumToAdd = 1
    LSeq = Len(StrainSeq(0))
    ReDim HitsPerseq(Nextno, LSeq)

    For X = 0 To Nextno
        EY = CurrentXOver(X)

        For Y = 1 To EY

            If XOverList(X, Y).Ending > LSeq Then
                XOverList(X, Y).Ending = LSeq
                XOverList(X, Y).Beginning = LSeq - 15
            End If

            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            BPos = XOverList(X, Y).Beginning
            EPos = XOverList(X, Y).Ending

            If BPos <= EPos Then
                SubM = Int((EPos - BPos) / 3)
                BPos = BPos + SubM
                EPos = EPos - SubM
                
                AddScores BPos, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
'                For Z = BPos To EPos
'                    HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
'                    'HitsPerSeq(Ma, Z) = HitsPerSeq(Ma, Z) + 1
'                    'HitsPerSeq(Mi, Z) = HitsPerSeq(Mi, Z) + 1
'                Next 'Z

            Else
                SubM = Int((EPos + (LSeq - BPos)) / 3)
                BPos = BPos + SubM

                If BPos > LSeq Then BPos = BPos - LSeq
                EPos = EPos - SubM

                If EPos < 0 Then EPos = LSeq + EPos

                If EPos < BPos Then
                    AddScores 1, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    AddScores BPos, LSeq, Nextno, NumToAdd, X, HitsPerseq(0, 0)
           '         For Z = 1 To EPos
           '             HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
           '         Next 'Z'

'                    For Z = BPos To LSeq
'                        HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
'                    Next 'Z

                Else
                    AddScores BPos, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    
              '      For Z = BPos To EPos
              '          HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
              '      Next 'Z

                End If

            End If

        Next 'Y

    Next 'X
    For X = 0 To Nextno
        Y = 1
        Do While Y <= CurrentXOver(X)
            ReDim NumHits(2)
            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            

            If Mi <= Nextno And Mi <> X Then
                
                BPos = XOverList(X, Y).Beginning
                EPos = XOverList(X, Y).Ending
                
               ' If BPos <= EPos Then
               '     SubM = Int((EPos - BPos) / 3)
               '     BPos = BPos + SubM
               '     EPos = EPos - SubM
               ' Else
               '     SubM = Int((EPos + (LSeq - BPos)) / 3)
               '     BPos = BPos + SubM
   '
     '               If BPos > LSeq Then BPos = BPos - LSeq
     '               EPos = EPos - SubM
   '
   '                 If EPos < 0 Then EPos = LSeq + EPos
   '             End If
                
                If BPos <= EPos Then

                    'First count the number of hits for this sequence
                    DoHits2 BPos, EPos, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
       '             For Z = BPos To EPos
       '                 NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
       '                 NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
       '                 NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
       '             Next 'Z

                Else
                    
                    
                    ' if region overlaps the ends
                    'First count the number of hits for this sequence
                    DoHits2 1, EPos, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                    DoHits2 BPos, LSeq, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                    
     '               For Z = BPos To LSeq
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

     '               For Z = 1 To EPos
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

                End If

                If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                    Call UpdateXOList(Mi, CurrentXOver(), XOverList())
                    
                    XOverList(Mi, CurrentXOver(Mi)) = XOverList(X, Y)
                    XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                    XOverList(Mi, CurrentXOver(Mi)).Daughter = XOverList(X, Y).MinorP
                    XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(X, Y).Daughter
                    XOverList(Mi, CurrentXOver(Mi)).MajorP = XOverList(X, Y).MajorP

                    If XOverList(X, Y).OutsideFlag = 0 Then
                        XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                    ElseIf XOverList(X, Y).OutsideFlag = 1 Then
                        XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                    ElseIf XOverList(X, Y).OutsideFlag = 2 Then
                        temp = XOverList(Mi, CurrentXOver(Mi)).MinorP
                        XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(ActiveMinorP, CurrentXOver(ActiveMinorP)).MajorP
                        XOverList(Mi, CurrentXOver(Mi)).MajorP = temp
                        XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 0
                        XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                    Else
                        XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = XOverList(X, Y).OutsideFlag
                    End If
                    
                    If Y = CurrentXOver(X) Then
                        CurrentXOver(X) = CurrentXOver(X) - 1
                    ElseIf CurrentXOver(X) > 0 Then
                        XOverList(X, Y) = XOverList(X, CurrentXOver(X))
                        CurrentXOver(X) = CurrentXOver(X) - 1
                        Y = Y - 1
                    End If

                ElseIf NumHits(2) > NumHits(0) And NumHits(2) > NumHits(1) Then
                    Call UpdateXOList(Ma, CurrentXOver(), XOverList())
                    XOverList(Ma, CurrentXOver(Ma)) = XOverList(X, Y)
                    XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                    XOverList(Ma, CurrentXOver(Ma)).Daughter = XOverList(X, Y).MajorP
                    XOverList(Ma, CurrentXOver(Ma)).MinorP = XOverList(X, Y).MinorP
                    XOverList(Ma, CurrentXOver(Ma)).MajorP = XOverList(X, Y).Daughter
                    XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2

                   
                    If XOverList(X, Y).OutsideFlag = 1 Then
                        XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                    ElseIf XOverList(X, Y).OutsideFlag = 2 Then
                        temp = XOverList(Ma, CurrentXOver(Ma)).MinorP
                        XOverList(Ma, CurrentXOver(Ma)).MinorP = XOverList(Ma, CurrentXOver(Ma)).MajorP
                        XOverList(Ma, CurrentXOver(Ma)).MajorP = temp
                        XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 0
                        XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                    ElseIf XOverList(X, Y).OutsideFlag = 0 Then
                        XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                        XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                    Else
                        XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = XOverList(X, Y).OutsideFlag
                    End If
                    
                    If Y = CurrentXOver(X) Then
                        CurrentXOver(X) = CurrentXOver(X) - 1
                    ElseIf CurrentXOver(X) > 1 Then
                        XOverList(X, Y) = XOverList(X, CurrentXOver(X))
                        CurrentXOver(X) = CurrentXOver(X) - 1
                        Y = Y - 1
                    
                    End If

                End If

            End If

            Y = Y + 1
        Loop

    Next 'X

End Sub
Public Sub ReCheck2()
    Dim HitsPerseq() As Long, HitsPerSeq2() As Long
    Dim NumHits() As Long
    Dim LSeq As Long, BPos As Long, EPos As Long, Ma As Long, Mi As Long, SubM As Long
    Dim X As Long, EY As Long, Y As Long, z As Long
    Dim NumToAdd As Long
    NumToAdd = 1
    LSeq = Len(StrainSeq(0))
    
    ReDim HitsPerseq(Nextno, LSeq)
    ReDim HitsPerSeq2(Nextno, LSeq)
    'Exit Sub
    For X = 0 To Nextno
        EY = CurrentXOver(X)

        For Y = 1 To EY

            If XOverList(X, Y).Ending > LSeq Then
                XOverList(X, Y).Ending = LSeq
                XOverList(X, Y).Beginning = LSeq - 15
            End If

            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            BPos = XOverList(X, Y).Beginning
            EPos = XOverList(X, Y).Ending

            If BPos <= EPos Then
                SubM = Int((EPos - BPos) / 3)
                BPos = BPos + SubM
                EPos = EPos - SubM
                
                AddScores BPos, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
'                For Z = BPos To EPos
'                    HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
'                    'HitsPerSeq(Ma, Z) = HitsPerSeq(Ma, Z) + 1
'                    'HitsPerSeq(Mi, Z) = HitsPerSeq(Mi, Z) + 1
'                Next 'Z

            Else
                SubM = Int((EPos + (LSeq - BPos)) / 3)
                BPos = BPos + SubM

                If BPos > LSeq Then BPos = BPos - LSeq
                EPos = EPos - SubM

                If EPos < 0 Then EPos = LSeq + EPos

                If EPos < BPos Then
                    AddScores 1, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    AddScores BPos, LSeq, Nextno, NumToAdd, X, HitsPerseq(0, 0)
           '         For Z = 1 To EPos
           '             HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
           '         Next 'Z'

'                    For Z = BPos To LSeq
'                        HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
'                    Next 'Z

                Else
                    AddScores BPos, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    
              '      For Z = BPos To EPos
              '          HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + 1
              '      Next 'Z

                End If

            End If

        Next 'Y
        Form1.ProgressBar1.Value = ((X / Nextno) * 0.2) * 100
    Next 'X
    
    For X = 0 To Nextno
        Y = 1
        Do While Y <= CurrentXOver(X)
            ReDim NumHits(2)
            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            

            If Mi <= Nextno And Mi <> X Then
                
                BPos = XOverList(X, Y).Beginning
                EPos = XOverList(X, Y).Ending
                
                If BPos <= EPos Then
                    SubM = Int((EPos - BPos) / 4)
                    BPos = BPos + SubM
                    EPos = EPos - SubM
                Else
                    SubM = Int((EPos + (LSeq - BPos)) / 4)
                    BPos = BPos + SubM
   
                   If BPos > LSeq Then BPos = BPos - LSeq
                    EPos = EPos - SubM
   
                    If EPos < 0 Then EPos = LSeq + EPos
                End If
                
                If BPos <= EPos Then

                    'First count the number of hits for this sequence
                    DoHits2 BPos, EPos, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
       '             For Z = BPos To EPos
       '                 NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
       '                 NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
       '                 NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
       '             Next 'Z
                    If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                        For z = BPos To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) + 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                        Next 'Z
                    ElseIf NumHits(2) > NumHits(0) And (NumHits(2) > NumHits(1)) Then
                        For z = BPos To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) + 1
                            
                        Next 'Z
                    ElseIf NumHits(0) > NumHits(1) And (NumHits(0) > NumHits(2)) Then
                        For z = BPos To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                        Next 'Z
                    Else
                      '  For Z = BPos To EPos
                      '      HitsPerSeq2(X, Z) = HitsPerSeq2(X, Z) - 1
                       '     HitsPerSeq2(Mi, Z) = HitsPerSeq2(Mi, Z) - 1
                      '      HitsPerSeq2(Ma, Z) = HitsPerSeq2(Ma, Z) - 1
                      '  Next 'Z
                    End If

                Else
                    
                    
                    ' if region overlaps the ends
                    'First count the number of hits for this sequence
                    DoHits 1, EPos, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                    DoHits BPos, LSeq, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                    
                    If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                        For z = 1 To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) + 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                        Next 'Z
                        For z = BPos To LSeq
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) + 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                        Next 'Z
                    ElseIf NumHits(2) > NumHits(0) And (NumHits(2) > NumHits(1)) Then
                        For z = 1 To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) + 1
                        Next 'Z
                        For z = BPos To LSeq
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) + 1
                        Next 'Z
                    ElseIf NumHits(0) > NumHits(1) And (NumHits(0) > NumHits(2)) Then
                        For z = 1 To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                        Next 'Z
                        For z = BPos To LSeq
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 1
                            HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                            HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                        Next 'Z
                    End If
                    
     '               For Z = BPos To LSeq
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

     '               For Z = 1 To EPos
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

                End If
                
            End If
            Y = Y + 1
        Loop
        Form1.ProgressBar1.Value = 20 + ((X / Nextno) * 0.2) * 100
    Next 'X
                
    For X = 0 To Nextno
        
        Y = 1
        Do While Y <= CurrentXOver(X)
            
          
            ReDim NumHits(2)
            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            

            If Mi <= Nextno And Mi <> X Then
                
                BPos = XOverList(X, Y).Beginning
                EPos = XOverList(X, Y).Ending
                
                If BPos <= EPos Then
                    SubM = Int((EPos - BPos) / 4)
                    BPos = BPos + SubM
                    EPos = EPos - SubM
                Else
                    SubM = Int((EPos + (LSeq - BPos)) / 4)
                    BPos = BPos + SubM
   
                   If BPos > LSeq Then BPos = BPos - LSeq
                    EPos = EPos - SubM
   
                    If EPos < 0 Then EPos = LSeq + EPos
                End If
                
                If BPos <= EPos Then

                    'First count the number of hits for this sequence
                    DoHits2 BPos, EPos, Nextno, X, Mi, Ma, HitsPerSeq2(0, 0), NumHits(0)


                Else
                    
                    
                    ' if region overlaps the ends
                    'First count the number of hits for this sequence
                    DoHits2 1, EPos, Nextno, X, Mi, Ma, HitsPerSeq2(0, 0), NumHits(0)
                    DoHits2 BPos, LSeq, Nextno, X, Mi, Ma, HitsPerSeq2(0, 0), NumHits(0)
     '               For Z = BPos To LSeq
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

     '               For Z = 1 To EPos
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

                End If
                If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                    Call UpdateXOList(Mi, CurrentXOver(), XOverList())
                    
                    XOverList(Mi, CurrentXOver(Mi)) = XOverList(X, Y)
                    'XOverList(Mi, CurrentXover(Mi)).MissIdentifyFlag = 1
                    XOverList(Mi, CurrentXOver(Mi)).Daughter = XOverList(X, Y).MinorP
                    XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(X, Y).Daughter
                    XOverList(Mi, CurrentXOver(Mi)).MajorP = XOverList(X, Y).MajorP
                    
                    d = X
                    Outer1 = 0
                    If XOverList(X, Y).MissIdentifyFlag = 1 Then
                    
                        If Distance(Mi, Ma) > (Distance(Mi, d) + Distance(Ma, d)) / 2 Then
                            Outer1 = 1
                        End If
                        
                    End If
                    If XOverList(X, Y).OutsideFlag = 0 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = XOverList(X, Y).OutsideFlag
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 2
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 1 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                            temp = XOverList(Mi, CurrentXOver(Mi)).MinorP
                            XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(Mi, CurrentXOver(Mi)).MajorP
                            XOverList(Mi, CurrentXOver(Mi)).MajorP = temp
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                            
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 2 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            Outer1 = 0
                            If Distance(Mi, Ma) < (Distance(Mi, d) + Distance(Ma, d)) / 2 Then
                                Outer1 = 1
                            End If
                            
                            If Outer1 = 1 And X = 1233456 Then
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                            Else
                                temp = XOverList(Mi, CurrentXOver(Mi)).MinorP
                                XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(Mi, CurrentXOver(Mi)).MajorP
                                XOverList(Mi, CurrentXOver(Mi)).MajorP = temp
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 1
                                XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                            End If
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            If Outer1 = 1 Or X = X Then
                                                                
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 1
                                
                            Else
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                            End If
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                        End If
                    End If
                    
                    If Y = CurrentXOver(X) Then
                        CurrentXOver(X) = CurrentXOver(X) - 1
                    ElseIf CurrentXOver(X) > 0 Then
                        XOverList(X, Y) = XOverList(X, CurrentXOver(X))
                        
                        CurrentXOver(X) = CurrentXOver(X) - 1
                        Y = Y - 1
                    End If

                ElseIf NumHits(2) > NumHits(0) And NumHits(2) > NumHits(1) Then
                    Call UpdateXOList(Ma, CurrentXOver(), XOverList())
                    
                    XOverList(Ma, CurrentXOver(Ma)).Beginning = XOverList(X, Y).Beginning
                    'XOverList(Ma, CurrentXover(Ma)).MissIdentifyFlag = 1
                    XOverList(Ma, CurrentXOver(Ma)).Daughter = XOverList(X, Y).MajorP
                    XOverList(Ma, CurrentXOver(Ma)).MinorP = XOverList(X, Y).MinorP
                    XOverList(Ma, CurrentXOver(Ma)).MajorP = XOverList(X, Y).Daughter
                    d = X
                    Outer1 = 0
                    If XOverList(X, Y).MissIdentifyFlag = 1 Then
                    
                        If Distance(Mi, Ma) > (Distance(Mi, d) + Distance(Ma, d)) / 2 Then
                            Outer1 = 1
                        End If
                        
                    End If
                    If XOverList(X, Y).OutsideFlag = 0 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = XOverList(X, Y).OutsideFlag
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 2
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 1 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            temp = XOverList(Ma, CurrentXOver(Ma)).MinorP
                            XOverList(Ma, CurrentXOver(Ma)).MinorP = XOverList(Ma, CurrentXOver(Ma)).MajorP
                            XOverList(Ma, CurrentXOver(Ma)).MajorP = temp
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                            
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 2 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 0
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            If Outer1 = 1 Then
                                                                
                                XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 1
                                
                            Else
                                XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            End If
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                        End If
                    End If
                    
                    If Y = CurrentXOver(X) Then
                        CurrentXOver(X) = CurrentXOver(X) - 1
                    ElseIf CurrentXOver(X) > 1 Then
                        XOverList(X, Y) = XOverList(X, CurrentXOver(X))
                        CurrentXOver(X) = CurrentXOver(X) - 1
                        Y = Y - 1
                    
                    End If

                End If

            End If

            Y = Y + 1
        Loop
        Form1.ProgressBar1.Value = 40 + ((X / Nextno) * 0.6) * 100
    Next 'X
    
End Sub

Public Sub ReSortB()

    Dim HitsPerseq() As Long, HitsPerSeq2() As Long
    Dim NumHits() As Long
    Dim LSeq As Long, BPos As Long, EPos As Long, Ma As Long, Mi As Long, SubM As Long
    Dim NumToAdd As Long, X As Long, EY As Long, Y As Long, z As Long
    
    LSeq = Len(StrainSeq(0))
    ReDim HitsPerseq(Nextno, LSeq)
    ReDim HitsPerSeq2(Nextno, LSeq)
    Form1.SSPanel1.Caption = "Re-Sorting Regions"
    For X = 0 To Nextno
        EY = CurrentXOver(X)

        For Y = 1 To EY
            If XOverList(X, Y).Ending > LSeq Then
                XOverList(X, Y).Ending = LSeq
                XOverList(X, Y).Beginning = LSeq - 15
            End If

            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            BPos = XOverList(X, Y).Beginning
            EPos = XOverList(X, Y).Ending
            If XOverList(X, Y).Accept = 0 Then
                NumToAdd = 1
            ElseIf XOverList(X, Y).Accept = 2 Then
                NumToAdd = 0
            ElseIf XOverList(X, Y).Accept = 1 Then
                NumToAdd = 10000
            End If
            
            If BPos <= EPos Then
                SubM = Int((EPos - BPos) / 3)
                BPos = BPos + SubM
                EPos = EPos - SubM
                
                'AddScores BPos, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                
                For z = BPos To EPos
                    HitsPerseq(X, z) = HitsPerseq(X, z) + NumToAdd
                    X = X
                Next 'Z

            Else
                SubM = Int((EPos + (LSeq - BPos)) / 3)
                BPos = BPos + SubM

                If BPos > LSeq Then BPos = BPos - LSeq
                EPos = EPos - SubM

                If EPos < 0 Then EPos = LSeq + EPos

                If EPos < BPos Then
                    
                    AddScores 1, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    AddScores BPos, LSeq, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    'For Z = 1 To EPos
                    '    HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + NumToAdd
                    'Next 'Z
                    
                    'For Z = BPos To LSeq
                    '    HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + NumToAdd
                    'Next 'Z

                Else
                    AddScores BPos, EPos, Nextno, NumToAdd, X, HitsPerseq(0, 0)
                    'For Z = BPos To EPos
                    '    HitsPerSeq(X, Z) = HitsPerSeq(X, Z) + NumToAdd
                   ' Next 'Z

                End If

            End If

        Next 'Y
        Form1.ProgressBar1 = (X / Nextno) * 30
    Next 'X

    For X = 0 To Nextno
        Y = 1
        Do While Y <= CurrentXOver(X)
            
            ReDim NumHits(2)
            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            

            If Mi <= Nextno And Mi <> X Then
                
                BPos = XOverList(X, Y).Beginning
                EPos = XOverList(X, Y).Ending
                
                If BPos <= EPos Then
                    SubM = Int((EPos - BPos) / 3)
                    BPos = BPos + SubM
                    EPos = EPos - SubM
                Else
                    SubM = Int((EPos + (LSeq - BPos)) / 3)
                    BPos = BPos + SubM
   
                   If BPos > LSeq Then BPos = BPos - LSeq
                    EPos = EPos - SubM
   
                    If EPos < 0 Then EPos = LSeq + EPos
                End If
                
                If BPos <= EPos Then
                    If XOverList(X, Y).Accept = 1 Then
                        For z = BPos To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 10000
                        Next z
                    Else
                        'First count the number of hits for this sequence
                        DoHits2 BPos, EPos, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                        'For Z = BPos To EPos
                        '    NumHits(0) = NumHits(0) + HitsPerseq(X, Z)
                        '    NumHits(1) = NumHits(1) + HitsPerseq(Mi, Z)
                        '    NumHits(2) = NumHits(2) + HitsPerseq(Ma, Z)
                        'Next 'Z
                        If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                            For z = BPos To EPos
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) + 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                            Next 'Z
                        ElseIf NumHits(2) > NumHits(0) And (NumHits(2) > NumHits(1)) Then
                            For z = BPos To EPos
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) + 1
                            Next 'Z
                        ElseIf NumHits(0) > NumHits(1) And (NumHits(0) > NumHits(2)) Then
                            For z = BPos To EPos
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                            Next 'Z
                        Else
                          '  For Z = BPos To EPos
                          '      HitsPerSeq2(X, Z) = HitsPerSeq2(X, Z) - 1
                           '     HitsPerSeq2(Mi, Z) = HitsPerSeq2(Mi, Z) - 1
                          '      HitsPerSeq2(Ma, Z) = HitsPerSeq2(Ma, Z) - 1
                          '  Next 'Z
                        End If
                    End If
                Else
                    
                    If XOverList(X, Y).Accept = 1 Then
                        For z = 1 To EPos
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 10000
                        Next z
                        For z = BPos To LSeq
                            HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 10000
                        Next z
                    Else
                        ' if region overlaps the ends
                        'First count the number of hits for this sequence
                        DoHits 1, EPos, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                        DoHits BPos, LSeq, Nextno, X, Mi, Ma, HitsPerseq(0, 0), NumHits(0)
                        
                        If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                            For z = 1 To EPos
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) + 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                            Next 'Z
                            For z = BPos To LSeq
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) + 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                            Next 'Z
                        ElseIf NumHits(2) > NumHits(0) And (NumHits(2) > NumHits(1)) Then
                            For z = 1 To EPos
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) + 1
                            Next 'Z
                            For z = BPos To LSeq
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) - 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) + 1
                            Next 'Z
                        ElseIf NumHits(0) > NumHits(1) And (NumHits(0) > NumHits(2)) Then
                            For z = 1 To EPos
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                            Next 'Z
                            For z = BPos To LSeq
                                HitsPerSeq2(X, z) = HitsPerSeq2(X, z) + 1
                                HitsPerSeq2(Mi, z) = HitsPerSeq2(Mi, z) - 1
                                HitsPerSeq2(Ma, z) = HitsPerSeq2(Ma, z) - 1
                            Next 'Z
                        End If
                    End If
     '               For Z = BPos To LSeq
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

     '               For Z = 1 To EPos
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

                End If
                
            End If
            Y = Y + 1
        Loop
        Form1.ProgressBar1 = 30 + (X / Nextno) * 30
    Next 'X
                
    For X = 0 To Nextno
        Y = 1
        Do While Y <= CurrentXOver(X)
            ReDim NumHits(2)
            Mi = XOverList(X, Y).MinorP
            Ma = XOverList(X, Y).MajorP
            

            If Mi <= Nextno And Mi <> X Then
                
                BPos = XOverList(X, Y).Beginning
                EPos = XOverList(X, Y).Ending
                
                If BPos <= EPos Then
                    SubM = Int((EPos - BPos) / 3)
                    BPos = BPos + SubM
                    EPos = EPos - SubM
                Else
                    SubM = Int((EPos + (LSeq - BPos)) / 3)
                    BPos = BPos + SubM
   
                   If BPos > LSeq Then BPos = BPos - LSeq
                    EPos = EPos - SubM
   
                    If EPos < 0 Then EPos = LSeq + EPos
                End If
                
                If BPos <= EPos Then

                    'First count the number of hits for this sequence
                    DoHits2 BPos, EPos, Nextno, X, Mi, Ma, HitsPerSeq2(0, 0), NumHits(0)


                Else
                    
                    
                    ' if region overlaps the ends
                    'First count the number of hits for this sequence
                    DoHits2 1, EPos, Nextno, X, Mi, Ma, HitsPerSeq2(0, 0), NumHits(0)
                    DoHits2 BPos, LSeq, Nextno, X, Mi, Ma, HitsPerSeq2(0, 0), NumHits(0)
     '               For Z = BPos To LSeq
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

     '               For Z = 1 To EPos
     '                   NumHits(0) = NumHits(0) + HitsPerSeq(X, Z)
     '                   NumHits(1) = NumHits(1) + HitsPerSeq(Mi, Z)
     '                   NumHits(2) = NumHits(2) + HitsPerSeq(Ma, Z)
     '               Next 'Z

                End If
                If NumHits(1) > NumHits(0) And (NumHits(1) > NumHits(2)) Then
                    Call UpdateXOList(Mi, CurrentXOver(), XOverList())
                    
                    XOverList(Mi, CurrentXOver(Mi)) = XOverList(X, Y)
                    'XOverList(Mi, CurrentXover(Mi)).MissIdentifyFlag = 1
                    XOverList(Mi, CurrentXOver(Mi)).Daughter = XOverList(X, Y).MinorP
                    XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(X, Y).Daughter
                    XOverList(Mi, CurrentXOver(Mi)).MajorP = XOverList(X, Y).MajorP
                    Outer1 = 0
                    d = X
                    If XOverList(X, Y).MissIdentifyFlag = 1 Then
                    
                        If Distance(Mi, Ma) > (Distance(Mi, d) + Distance(Ma, d)) / 2 Then
                            Outer1 = 1
                        End If
                        
                    End If
                    If XOverList(X, Y).OutsideFlag = 0 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 1
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = XOverList(X, Y).OutsideFlag
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 1
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 2
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 1 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 0
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = XOverList(X, Y).OutsideFlag
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 2 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            temp = XOverList(Mi, CurrentXOver(Mi)).MinorP
                            XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(Mi, CurrentXOver(Mi)).MajorP
                            XOverList(Mi, CurrentXOver(Mi)).MajorP = temp
                            XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 0
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            If Outer1 = 1 Or X = X Then
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 1
                                temp = XOverList(Mi, CurrentXOver(Mi)).MinorP
                                XOverList(Mi, CurrentXOver(Mi)).MinorP = XOverList(Mi, CurrentXOver(Mi)).MajorP
                                XOverList(Mi, CurrentXOver(Mi)).MajorP = temp
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 0
                                
                            Else
                                XOverList(Mi, CurrentXOver(Mi)).OutsideFlag = 2
                            End If
                            XOverList(Mi, CurrentXOver(Mi)).MissIdentifyFlag = 1
                        End If
                    End If
                    
                    Call CompressList(X, Y, XOverList(), CurrentXOver())
                    
                    If CurrentXOver(X) > 1 And Y < CurrentXOver(X) Then Y = Y - 1
                ElseIf NumHits(2) > NumHits(0) And NumHits(2) > NumHits(1) Then
                    Call UpdateXOList(Ma, CurrentXOver(), XOverList())
                    XOverList(Ma, CurrentXOver(Ma)) = XOverList(X, Y)

                    'XOverList(Ma, CurrentXover(Ma)).MissIdentifyFlag = 1
                    XOverList(Ma, CurrentXOver(Ma)).Daughter = XOverList(X, Y).MajorP
                    XOverList(Ma, CurrentXOver(Ma)).MinorP = XOverList(X, Y).MinorP
                    XOverList(Ma, CurrentXOver(Ma)).MajorP = XOverList(X, Y).Daughter
                    d = X
                    Outer1 = 0
                    If XOverList(X, Y).MissIdentifyFlag = 1 Then
                    
                        If Distance(Mi, Ma) > (Distance(Mi, d) + Distance(Ma, d)) / 2 Then
                            Outer1 = 1
                        End If
                        
                    End If
                    If XOverList(X, Y).OutsideFlag = 0 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = XOverList(X, Y).OutsideFlag
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 2
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 1 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            temp = XOverList(Ma, CurrentXOver(Ma)).MinorP
                            XOverList(Ma, CurrentXOver(Ma)).MinorP = XOverList(Ma, CurrentXOver(Ma)).MajorP
                            XOverList(Ma, CurrentXOver(Ma)).MajorP = temp
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                            
                        End If
                    ElseIf XOverList(X, Y).OutsideFlag = 2 Then
                        If XOverList(X, Y).MissIdentifyFlag = 0 Or XOverList(X, Y).MissIdentifyFlag = 2 Then
                            
                            XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 0
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = XOverList(X, Y).MissIdentifyFlag
                        ElseIf XOverList(X, Y).MissIdentifyFlag = 1 Then
                            If Outer1 = 1 Then
                                                                
                                XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 1
                                
                            Else
                                XOverList(Ma, CurrentXOver(Ma)).OutsideFlag = 2
                            End If
                            XOverList(Ma, CurrentXOver(Ma)).MissIdentifyFlag = 1
                        End If
                    End If
                    
                    Call CompressList(X, Y, XOverList(), CurrentXOver())
                    
                    If CurrentXOver(X) > 1 And Y < CurrentXOver(X) Then
                        Y = Y - 1
                    
                    End If

                End If

            End If

            Y = Y + 1
        Loop
        Form1.ProgressBar1 = 60 + (X / Nextno) * 30
    Next 'X

End Sub

Public Sub LoadLF(FN As String, DoSeqNo)
Dim FNum As Integer
ErrorFlag = 0
FNum = FreeFile
Open "infile" For Output As #FNum
Print #FNum, FN
Print #FNum, Trim(Str(DoSeqNo + 1))
Close #FNum
ExpectFL = DoSeqNo
'OP = GetCommandOutput("lkgen.bat", 5, True, True)
Call ShellAndClose("lkgen.bat", 0)

If FN = "LF100" Then
   FN = "LF" + Trim(Str(DoSeqNo + 1))
End If
On Error Resume Next

Kill FN
On Error GoTo 0
Open "new_lk.txt" For Append As #FNum
If LOF(FNum) = 0 Then
    Response = MsgBox("There was a problem extracting data from the likelihood lookup file.  Make sure that files named 'lkgen.exe' and 'LF100' are in the RDP startup directory", 0, "RDP Warning")
    
    Close #FNum
    VRFlag = 0
    Form1.SSPanel1.Caption = ""
    ErrorFlag = 1
    Exit Sub
End If
Close #FNum
On Error Resume Next
    Kill "LF" & (Trim(Str(DoSeqNo + 1)))
On Error GoTo 0

Name "new_lk.txt" As "LF" & (Trim(Str(DoSeqNo + 1)))

FN = "LF" + Trim(Str(DoSeqNo + 1))
End Sub
Public Sub MakePVO5(RS, PVO As String)
Pos = InStr(1, PVO, ".", vbBinaryCompare) 'Check for decimal

If Pos = 0 Then
    If Len(PVO) < RS - 1 Then
        PVO = PVO + "."
    Else
        Exit Sub
    End If
End If
If Len(PVO) = RS Then
    Exit Sub
ElseIf Len(PVO) = RS - 1 Then
    PVO = PVO + "0"
ElseIf Len(PVO) = RS - 2 Then
    PVO = PVO + "00"
ElseIf Len(PVO) = RS - 3 Then
    PVO = PVO + "000"
ElseIf Len(PVO) = RS - 4 Then
    PVO = PVO + "0000"
ElseIf Len(PVO) = RS - 5 Then
    PVO = PVO + "00000"
ElseIf Len(PVO) = RS - 6 Then
    PVO = PVO + "00000"
ElseIf Len(PVO) = RS - 7 Then
    PVO = PVO + "00000"
End If

End Sub
Public Sub NextSeqs(GoOn)
GoOn = 1
Seq3 = Seq3 + 1
If Seq3 > Nextno Then
    Seq2 = Seq2 + 1
    If Seq2 < Nextno Then
        Seq3 = Seq2 + 1
    Else
        Seq1 = Seq1 + 1
        If Seq1 < Nextno - 1 Then
            Seq2 = Seq1 + 1
            Seq3 = Seq2 + 1
        Else
            GoOn = 0
        End If
    End If
End If
If Seq1 > Nextno Or Seq2 > Nextno Or Seq3 > Nextno Then
    X = X
End If

End Sub
Public Sub RefreshTimes()
Dim oShowPlt As Byte
If DontRefreshFlag = 1 Then Exit Sub
oShowPlt = ShowPlotFlag
ShowPlotFlag = 0
If NoF3Check = 1 Or NoF3Check2 = 1 Then
    ShowPlotFlag = oShowPlt
    Exit Sub
End If
If Form3.Visible = True And Form3.Command1.Enabled = True And Form3.Command1.Visible = True Then
    Form3.Command1.SetFocus
End If
Dim TotT As Double
    
If Form3.Check1.Value = 1 Then
    DoScans(0, 2) = 1
Else
    DoScans(0, 2) = 0
End If
If Form3.Check12.Value = 1 Then
    DoScans(0, 8) = 1
Else
    DoScans(0, 8) = 0
End If
If Form3.Check2.Value = 1 Then
    DoScans(0, 3) = 1
Else
    DoScans(0, 3) = 0
End If
If Form3.Check3.Value = 1 Then
    DoScans(0, 4) = 1
Else
    DoScans(0, 4) = 0
End If
If Form3.Check4.Value = 1 Then
    DoScans(0, 0) = 1
Else
    DoScans(0, 0) = 0
End If
If Form3.Check5.Value = 1 Then
    DoScans(0, 1) = 1
Else
    DoScans(0, 1) = 0
End If
If Form3.Check6.Value = 1 Then
    DoScans(0, 5) = 1
Else
    DoScans(0, 5) = 0
End If
If Form3.Check12.Value = 1 Then
    DoScans(0, 8) = 1
Else
    DoScans(0, 8) = 0
End If
If Form3.Check21.Value = 1 Then
    DoScans(0, 7) = 1
Else
    DoScans(0, 7) = 0
End If



If Form3.Check14.Value = 1 Then
    DoScans(1, 0) = 1
Else
    DoScans(1, 0) = 0
End If
If Form3.Check15.Value = 1 Then
    DoScans(1, 1) = 1
Else
    DoScans(1, 1) = 0
End If
If Form3.Check16.Value = 1 Then
    DoScans(1, 2) = 1
Else
    DoScans(1, 2) = 0
End If
If Form3.Check17.Value = 1 Then
    DoScans(1, 3) = 1
Else
    DoScans(1, 3) = 0
End If
If Form3.Check18.Value = 1 Then
    DoScans(1, 4) = 1
Else
    DoScans(1, 4) = 0
End If
If Form3.Check19.Value = 1 Then
    DoScans(1, 5) = 1
Else
    DoScans(1, 5) = 0
End If
If Form3.Check20.Value = 1 Then
    DoScans(1, 8) = 1
Else
    DoScans(1, 8) = 0
End If
If Form3.Check22.Value = 1 Then
    DoScans(1, 7) = 1
Else
    DoScans(1, 7) = 0
End If

Call GetTot(TotT)

Call DoProgLine(TotT)

ShowPlotFlag = oShowPlt
End Sub

Public Sub SetChecks()
NoF3Check = 1
If DoScans(0, 2) = 1 Then
    Form3.Check1.Value = 1
Else
    Form3.Check1.Value = 0
End If

If DoScans(0, 8) = 1 Then
    Form3.Check12.Value = 1
Else
    Form3.Check12.Value = 0
End If

If DoScans(0, 3) = 1 Then
    Form3.Check2.Value = 1
Else
    Form3.Check2.Value = 0
End If

If DoScans(0, 4) = 1 Then
    Form3.Check3.Value = 1
Else
    Form3.Check3.Value = 0
End If

If DoScans(0, 0) = 1 Then
    Form3.Check4.Value = 1
Else
    Form3.Check4.Value = 0
End If

If DoScans(0, 1) = 1 Then
    Form3.Check5.Value = 1
Else
    Form3.Check5.Value = 0
End If

If DoScans(0, 5) = 1 Then
    Form3.Check6.Value = 1
Else
    Form3.Check6.Value = 0
End If

If DoScans(0, 8) = 1 Then
    Form3.Check12.Value = 1
Else
    Form3.Check12.Value = 0
End If

If DoScans(0, 7) = 1 Then
    Form3.Check21.Value = 1
Else
    Form3.Check21.Value = 0
End If

If DoScans(1, 0) = 1 Then
    Form3.Check14.Value = 1
Else
    Form3.Check14.Value = 0
End If

If DoScans(1, 1) = 1 Then
    Form3.Check15.Value = 1
Else
    Form3.Check15.Value = 0
End If

If DoScans(1, 2) = 1 Then
    Form3.Check16.Value = 1
Else
    Form3.Check16.Value = 0
End If

If DoScans(1, 3) = 1 Then
    Form3.Check17.Value = 1
Else
    Form3.Check17.Value = 0
End If

If DoScans(1, 4) = 1 Then
    Form3.Check18.Value = 1
Else
    Form3.Check18.Value = 0
End If

If DoScans(1, 5) = 1 Then
    Form3.Check19.Value = 1
Else
    Form3.Check19.Value = 0
End If
If DoScans(1, 8) = 1 Then
    Form3.Check20.Value = 1
Else
    Form3.Check20.Value = 0
End If
If DoScans(1, 7) = 1 Then
    Form3.Check22.Value = 1
Else
    Form3.Check22.Value = 0
End If
NoF3Check = 0
End Sub
Public Sub UpdateAgeScore(z, CurAge, Curscore, Beginning, Ending, AgeScore() As Double, EventScore() As Long)
    
Dim A As Long
    If Beginning < Ending Then
        For A = Beginning To Ending
            
            If AgeScore(A, z) > CurAge Or AgeScore(A, z) = -1 Then
                AgeScore(A, z) = CurAge
                EventScore(A, z) = Curscore
            End If
        Next A
    Else
        For A = 1 To Ending
            If AgeScore(A, z) > CurAge Or AgeScore(A, z) = -1 Then
                AgeScore(A, z) = CurAge
                EventScore(A, z) = Curscore
            End If
        Next A
        For A = Beginning To Len(StrainSeq(0))
            If AgeScore(A, z) > CurAge Or AgeScore(A, z) = -1 Then
                AgeScore(A, z) = CurAge
                EventScore(A, z) = Curscore
            End If
        Next A
    End If
End Sub
Public Sub SetUpScanArrays()

'If DoScans(0, 0) = 1 Then
    ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
    ReDim XDiffpos(Len(StrainSeq(0)) + 200)
    ReDim XPosDiff(Len(StrainSeq(0)) + 200)
    ReDim ValidSpacer(Nextno)
    ReDim SpacerSeqs(Nextno)
    ReDim XOverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)
'End If
    
'If DoScans(0, 1) = 1 Then
    'If GCtripletflag = 1 Then
    ReDim SubSeq(Len(StrainSeq(0)), 6)
    ReDim XDiffpos(Len(StrainSeq(0)) + 200)
    ReDim XPosDiff(Len(StrainSeq(0)) + 200)
    ReDim FragMaxScore(Len(StrainSeq(0)), 5)
    ReDim MaxScorePos(Len(StrainSeq(0)), 5)
    ReDim PVals(Len(StrainSeq(0)), 5)
    ReDim FragSt(Len(StrainSeq(0)), 6)
    ReDim FragEn(Len(StrainSeq(0)), 6)
    ReDim FragScore(Len(StrainSeq(0)), 6)
   ' End If
'End If
    
'If DoScans(0, 3) = 1 Then
    'If MCTripletFlag = 0 Then
    HWindowWidth = CLng(MCWinSize / 2)
    ReDim Scores(Len(StrainSeq(0)), 2)  ' 0=s1,s2Matches etc
    ReDim Winscores(Len(StrainSeq(0)) + HWindowWidth * 2, 2) ' 0=s1,s2Matches etc
    ReDim ChiVals(Len(StrainSeq(0)), 2)
    ReDim ChiPvals(Len(StrainSeq(0)), 2)
    ReDim SmoothChi(Len(StrainSeq(0)), 2)
    ReDim XDiffpos(Len(StrainSeq(0)) + 200)
    ReDim XPosDiff(Len(StrainSeq(0)) + 200)
    If MCProportionFlag = 0 Then
        
        Call GetCriticalDiff(0)
        
        If MCWinSize <> HWindowWidth * 2 And MCProportionFlag = 0 Then
            MCWinSize = HWindowWidth * 2
        End If
    End If
   ' End If
'End If
'If DoScans(0, 4) = 1 Then
    HWindowWidth = CLng(CWinSize / 2)
    ReDim ScoresX(Len(StrainSeq(0)))  ' 0=s1,s2Matches etc
    ReDim WinScoresX(Len(StrainSeq(0)) + HWindowWidth * 2) ' 0=s1,s2Matches etc
    ReDim ChiValsX(Len(StrainSeq(0)))
    ReDim ChiPValsX(Len(StrainSeq(0)))
    ReDim SmoothChiX(Len(StrainSeq(0)))
    ReDim XDiffpos(Len(StrainSeq(0)) + 200)
    ReDim XPosDiff(Len(StrainSeq(0)) + 200)
    Call GetCriticalDiff(1)
    If CWinSize <> HWindowWidth * 2 And CProportionFlag = 0 Then
        CWinSize = HWindowWidth * 2
    End If
    
'End If
    
    If XTableFlag = 0 Then
        Call Build3SeqTable
    End If
    ReDim XoverSeqNumTS(Len(StrainSeq(0)))

End Sub

Public Sub ProcessEvent(PRN, ProbabilityXOver As Double, BTarget As Long, ETarget As Long, XPosDiff() As Long, XDiffpos() As Long, BWarn, EWarn, XOverWindow, LenXoverSeq)

'AD = Daughter
If CurrentXOver(Seq1) <= CurrentXOver(Seq2) And CurrentXOver(Seq1) <= CurrentXOver(Seq3) Then
    AD = Seq1
    AMi = Seq2
    AMa = Seq3
ElseIf CurrentXOver(Seq2) <= CurrentXOver(Seq1) And CurrentXOver(Seq2) <= CurrentXOver(Seq3) Then
    AD = Seq2
    AMi = Seq1
    AMa = Seq3
Else
    AD = Seq3
    AMi = Seq1
    AMa = Seq2
End If
'Keep track of signal numbers
oRecombNo(100) = oRecombNo(100) + 1
oRecombNo(PRN) = oRecombNo(PRN) + 1
If APermFlag = 0 Then
    Call UpdateXOList3(AD, CurrentXOver(), XOverList(), PRN, ProbabilityXOver, SIP)
Else
    SIP = 1
End If
                                        
                                            
If MCFlag = 2 Then
    If -Log10(ProbabilityXOver) * 2 > 0 And -Log10(ProbabilityXOver) * 2 < 100 Then
        PValCat(CurrentCorrect, CInt(-Log10(ProbabilityXOver) * 2)) = PValCat(CurrentCorrect, CInt(-Log10(ProbabilityXOver) * 2)) + 1
    ElseIf CInt(-Log10(ProbabilityXOver) * 2) >= 100 Then
        PValCat(CurrentCorrect, 100) = PValCat(CurrentCorrect, 100) + 1
    End If
End If

'if p high eough then add it to list, if not discard then repeat from Z to en.
If SIP > 0 Then
    Call CentreBP(BTarget, ETarget, XPosDiff(), XDiffpos(), BWarn, EWarn, XOverWindow, LenXoverSeq)
    XOverList(AD, SIP).Beginning = BTarget
    XOverList(AD, SIP).Ending = ETarget
    
    XOverList(AD, SIP).MajorP = AMa
    XOverList(AD, SIP).MinorP = AMi
    XOverList(AD, SIP).Daughter = AD
    XOverList(AD, SIP).ProgramFlag = PRN
    XOverList(AD, SIP).Probability = ProbabilityXOver
    
    
    If SEventNumber = 0 And ShowPlotFlag = 2 And (CLine = "" Or CLine = " ") Then
         StartPlt(PRN) = 1
         Call UpdatePlotB(AD, AMa, AMi, SIP)
    
    End If
    'Make a reminder that one of the berakpoints went undetected
    If SEventNumber > 0 Then
        If EWarn = 0 Then Call CheckEndsVB(XOverWindow, EWarn, LenXoverSeq, 1, CircularFlag, Seq1, Seq2, Seq3, BTarget, ETarget, SeqNum(), XPosDiff(), XDiffpos())
        If BWarn = 0 Then Call CheckEndsVB(XOverWindow, BWarn, LenXoverSeq, 0, CircularFlag, Seq1, Seq2, Seq3, BTarget, ETarget, SeqNum(), XPosDiff(), XDiffpos())
    End If
    If BWarn = 1 And EWarn = 1 Then
        XOverList(AD, SIP).SBPFlag = 3
    ElseIf BWarn = 1 Then
        XOverList(AD, SIP).SBPFlag = 1
    ElseIf EWarn = 1 Then
        XOverList(AD, SIP).SBPFlag = 2
    End If

    
    
ElseIf DoneRedo = 0 Then
    DoneRedo = 1
    Call AddToRedoList(PRN, Seq1, Seq2, Seq3)
End If



End Sub

Public Sub GetFromFile(PValCon, AddNumX, SNNextNo, FN, X, Y, XOverList() As XOverDefine)
Dim TS As String

Call GetTS(FN, TS)
XOverList(X, Y).Daughter = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).MajorP = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).MinorP = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).Beginning = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).Ending = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).Probability = CDbl(TS)
Call GetTS(FN, TS)
XOverList(X, Y).OutsideFlag = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).MissIdentifyFlag = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).PermPVal = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).ProgramFlag = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).LHolder = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).DHolder = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).BeginP = CDbl(TS)
Call GetTS(FN, TS)
XOverList(X, Y).EndP = CDbl(TS)
Call GetTS(FN, TS)
XOverList(X, Y).SBPFlag = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).Accept = CLng(TS)
Call GetTS(FN, TS)
XOverList(X, Y).Eventnumber = CLng(TS)
                
If XOverList(X, Y).ProgramFlag >= AddNumX Then XOverList(X, Y).Accept = 2
XOverList(X, Y).MajorP = XOverList(X, Y).MajorP + SNNextNo
XOverList(X, Y).MinorP = XOverList(X, Y).MinorP + SNNextNo
XOverList(X, Y).Daughter = XOverList(X, Y).Daughter + SNNextNo
XOverList(X, Y).Probability = XOverList(X, Y).Probability * PValCon
If XOverList(X, Y).DHolder >= 0 Then
    XOverList(X, Y).DHolder = XOverList(X, Y).DHolder + SNNextNo
Else
    XOverList(X, Y).DHolder = (Abs(XOverList(X, Y).DHolder) + SNNextNo) * -1
End If
If XOverList(X, Y).ProgramFlag >= AddNumX Then
    XOverList(X, Y).ProgramFlag = XOverList(X, Y).ProgramFlag - AddNumX + AddNum
End If

End Sub
Public Sub GetTS(FN, TS As String)
On Error Resume Next
TS = ""
Do While Not EOF(FN)
    Input #FN, TS
    If TS <> "" Then Exit Sub
Loop
On Error GoTo 0
End Sub
Public Sub DrawLRDMat()

Form1.Picture26.ScaleMode = 3
Form1.Picture26.AutoRedraw = True
CurMatrixFlag = 11
DoneMatX(11) = 1
MatFlag = 11
MaxN = MatBound(MatFlag)
RSize = UBound(MatrixL, 1) - 1

Dim PosS(1) As Long, PosE(1) As Long, DistD As Long

XAddj = (Form1.Picture26.ScaleHeight) / RSize
DistD = RSize / MatZoom(MatFlag)

PosS(0) = MatCoord(MatFlag, 0)
PosE(0) = PosS(0) + DistD
PosS(1) = MatCoord(MatFlag, 1)
PosE(1) = PosS(1) + DistD

Dim Limit As Long
Limit = RSize + 1

    
If PosE(1) > (Limit - 1) - 1 Then PosE(1) = (Limit - 1) - 1
If PosE(0) > (Limit - 1) - 1 Then PosE(0) = (Limit - 1) - 1
       
If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
        

Form1.Picture26.Picture = LoadPicture()
Form1.Picture26.ScaleMode = 3

DistD = RSize / MatZoom(MatFlag)
XAddj = (Form1.Picture26.ScaleHeight) / DistD
     

Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixL(), HeatMap(), CurScale, MaxN)

If DontDoKey = 0 Then
    Call DoKey(1, MaxN, 0, MatFlag, "Log Likelihood Ratio", CurScale)
End If

BPos = XOverList(RelX, RelY).Beginning
EPos = XOverList(RelX, RelY).Ending


'BPCoord(0) = ((BPos - LRDStep / 2) / Len(StrainSeq(0))) * RSize
'BPCoord(1) = ((EPos + LRDStep / 2) / Len(StrainSeq(0))) * RSize
BPCoord(0) = ((BPos) / Len(StrainSeq(0))) * RSize
BPCoord(1) = ((EPos) / Len(StrainSeq(0))) * RSize

FirstPass = 1
Call DoSpot(11)
Form1.Picture26.Refresh
Form1.Picture17.Refresh
Form1.ProgressBar1.Value = 0
Form1.Picture26.AutoRedraw = tr
End Sub


Public Sub DoEnds(SPF)
Dim RL As Long, LS As Long, BPos As Long, EPos As Long
    If XOverList(RelX, RelY).ProgramFlag = 1 Then
            X = X
        End If
    If XOverList(RelX, RelY).ProgramFlag = 0 Then
        RL = XOverWindow
        
    ElseIf XOverList(RelX, RelY).ProgramFlag = 1 Then
        RL = 10
    ElseIf XOverList(RelX, RelY).ProgramFlag = 2 Then
        RL = 10
    ElseIf XOverList(RelX, RelY).ProgramFlag = 3 Then
        RL = MCWinSize / 2
    ElseIf XOverList(RelX, RelY).ProgramFlag = 4 Then
        RL = CWinSize / 2
    ElseIf XOverList(RelX, RelY).ProgramFlag = 5 Then
        RL = 10
    ElseIf XOverList(RelX, RelY).ProgramFlag = 8 Then
        RL = 10
    End If
    RL = 10
    LS = LenXoverSeq
    BPos = XOverList(RelX, RelY).Beginning
    EPos = XOverList(RelX, RelY).Ending
    If SPF = 0 Or SPF = 2 Then
        warn = 0
        Call CheckEndsVB(RL, warn, LS, 0, CircularFlag, Seq1, Seq2, Seq3, BPos, EPos, SeqNum(), XPosDiff(), XDiffpos())
        
        If warn = 1 Then
            If XOverList(RelX, RelY).SBPFlag = 0 Then
                XOverList(RelX, RelY).SBPFlag = 1
            ElseIf XOverList(RelX, RelY).SBPFlag = 2 Then
                XOverList(RelX, RelY).SBPFlag = 3
            End If
        Else
            If XOverList(RelX, RelY).SBPFlag = 1 Then
                XOverList(RelX, RelY).SBPFlag = 0
            ElseIf XOverList(RelX, RelY).SBPFlag = 3 Then
                XOverList(RelX, RelY).SBPFlag = 2
            End If
        End If
    End If
    If SPF = 1 Or SPF = 2 Then
        
        warn = 0
        Call CheckEndsVB(RL, warn, LS, 1, CircularFlag, Seq1, Seq2, Seq3, BPos, EPos, SeqNum(), XPosDiff(), XDiffpos())
    
        If warn = 1 Then
            If XOverList(RelX, RelY).SBPFlag = 0 Then
                XOverList(RelX, RelY).SBPFlag = 2
            ElseIf XOverList(RelX, RelY).SBPFlag = 1 Then
                XOverList(RelX, RelY).SBPFlag = 3
            End If
        Else
            If XOverList(RelX, RelY).SBPFlag = 2 Then
                XOverList(RelX, RelY).SBPFlag = 0
            ElseIf XOverList(RelX, RelY).SBPFlag = 3 Then
                XOverList(RelX, RelY).SBPFlag = 1
            End If
        End If
    End If
End Sub
Public Sub MCXoverK(SPF)



Dim TWin As Long, X As Long, LO As Long, RO As Long, N As Long, B As Long, d As Long, SPos As Long, EPos As Long, Step As Long, LSeq As Long, Last As Long, numSites As Long, NumPairs As Long, WindowWidth As Long, NumBreakpoints As Long
Dim lPrb As Double, mPrb As Double, MChi As Double, oMCWinSize As Long, oMCWinfract As Double

MatFlag = 8


LenXoverSeq = 0

SS = GetTickCount
If SPF = 1 Then
    Form1.ProgressBar1 = 5
End If
If DoneMatX(8) = 0 Then
    ReDim XDiffpos(Len(StrainSeq(0)) + 200)
    ReDim XPosDiff(Len(StrainSeq(0)) + 200)
    ReDim VarsitesMC(Len(StrainSeq(0)) + 200)
    LenXoverSeq = FindSubSeqC(Len(StrainSeq(0)) + 1, Nextno, Seq1, Seq2, Seq3, SeqNum(0, 0), VarsitesMC(0), XPosDiff(0))
    DoneMatX(8) = 1
    Form1.SSPanel1.Caption = "Drawing MaxChi breakpoint matrix"
    ReDim MatrixMC(LenXoverSeq + 1, LenXoverSeq + 1)
    If LenXoverSeq < 2 Then Exit Sub
        
    ReDim Scores(Len(StrainSeq(0)), 2)  ' 0=s1,s2Matches etc
    ReDim Winscores(Len(StrainSeq(0)) + HWindowWidth * 2, 2)
    
    Dim CurMChi As Double
    
    
    'ReDim LenFrag(LenXOverSeq, LenXOverSeq)
    'ReDim WinScores(LenXOverSeq + 1, LenXOverSeq + 1, 2, 1)
    
    HWindowWidth = 1
    Dummy = WinScoreCalc(CriticalDiff, HWindowWidth, LenXoverSeq, Len(StrainSeq(0)) + 1, Seq1, Seq2, Seq3, Scores(0, 0), VarsitesMC(0), SeqNum(0, 0), Winscores(0, 0))
    Dim A(2) As Double, C(2) As Double, CntX As Long, TargetCntX As Long
    If X = X Then
        Dummy = MakeMatrixMC(Len(StrainSeq(0)), LenXoverSeq, A(0), C(0), Scores(0, 0), MatrixMC(0, 0))
    Else
    
    
        CntX = 0: TargetCntX = LenXoverSeq * (LenXoverSeq - 1) / 2
        eee = 0
        For X = 1 To LenXoverSeq - 3
            For Y = X + 3 To LenXoverSeq
                CntX = CntX + 1
                If Y = X + 3 Then
                   CurMChi = 0
                   For z = 0 To 2
                       A(z) = 0
                       C(z) = 0
                       For P = 1 To X
                           A(z) = A(z) + Scores(P, z)
                       Next P
                       For P = X + 1 To Y
                           C(z) = C(z) + Scores(P, z)
                       Next P
                       For P = Y + 1 To LenXoverSeq
                           A(z) = A(z) + Scores(P, z)
                       Next P
                       B = (LenXoverSeq - (Y - X)) - A(z)
                        d = (Y - X) - (z)
                       thMChi = CalcChiV(A(z), B, C(z), d)
                       
                       If thMChi > CurMChi Then CurMChi = thMChi
                   Next z
                Else
                    CurMChi = 0
                    For z = 0 To 2
                        A(z) = A(z) - Scores(Y, z) '96
                        C(z) = C(z) + Scores(Y, z) '0
                        B = (LenXoverSeq - (Y - X)) - A(z) '757
                        d = (Y - X) - C(z) '3
                        thMChi = CalcChiV(A(z), B, C(z), d)
                        If thMChi > CurMChi Then CurMChi = thMChi
                    Next z
                End If
                
                MatrixMC(X, Y) = -Log10(ChiPVal(CurMChi))
                MatrixMC(Y, X) = MatrixMC(X, Y)
            Next Y
            SSS = GetTickCount
            If Abs(SSS - eee) > 500 Then
                eee = SSS
                Form1.ProgressBar1.Value = (CntX / TargetCntX) * 100
            
            End If
        Next X
    End If
    
    MaxN = FindMaxN(LenXoverSeq, MatrixMC(0, 0))
    'For X = 1 To LenXOverSeq
    '    For Y = X + 1 To LenXOverSeq
    '        If MatrixMC(X, Y) > MaxN Then MaxN = MatrixMC(X, Y)
    '    Next Y
    'Next X
    MatBound(8) = MaxN
Else
    MatFlag = 8
    MaxN = MatBound(MatFlag)
    LenXoverSeq = UBound(MatrixMC, 1) - 1
End If

If SPF = 1 Then
    Form1.ProgressBar1 = 10
End If

RSize = LenXoverSeq
Dim PosS(1) As Long, PosE(1) As Long, DistD As Long
XAddj = (Form1.Picture26.ScaleHeight) / RSize
DistD = RSize / MatZoom(MatFlag)
PosS(0) = MatCoord(MatFlag, 0)
PosE(0) = PosS(0) + DistD
PosS(1) = MatCoord(MatFlag, 1)
PosE(1) = PosS(1) + DistD
Dim Limit As Long
Limit = UBound(MatrixMC, 1)

    
If PosE(1) > (Limit - 1) - 1 Then PosE(1) = (Limit - 1) - 1
If PosE(0) > (Limit - 1) - 1 Then PosE(0) = (Limit - 1) - 1
       
If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
        

Form1.Picture26.Picture = LoadPicture()
Form1.Picture26.ScaleMode = 3
DistD = RSize / MatZoom(MatFlag)
XAddj = (Form1.Picture26.ScaleHeight) / DistD
 XX = Form1.Picture26.Top
SS = GetTickCount
Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixMC(), HeatMap(), CurScale, MaxN)
EE = GetTickCount
TT = EE - SS
If SPF = 1 Then
    Form1.ProgressBar1 = 95
End If
If DontDoKey = 0 Then
    Call DoKey(1, MaxN, 0, MatFlag, "-Log(Chi P-Val)", CurScale)
End If

BPos = XOverList(RelX, RelY).Beginning
EPos = XOverList(RelX, RelY).Ending
For X = 1 To LenXoverSeq - 1
    If VarsitesMC(X) <= BPos And VarsitesMC(X + 1) >= BPos Then
        BPCoord(0) = X
    End If
    If VarsitesMC(X) <= EPos And VarsitesMC(X + 1) >= EPos Then
        BPCoord(1) = X
    End If
    X = X
Next X
If SPF = 1 Then
    Form1.ProgressBar1 = 100
End If
FirstPass = 1
Call DoSpot(8)
Form1.Picture26.Refresh
Form1.Picture17.Refresh
Form1.ProgressBar1.Value = 0
Form1.Picture26.AutoRedraw = True
Form1.SSPanel1.Caption = ""
End Sub

Public Sub DoSpot(CurMatrixFlag)

If CurMatrixFlag = 8 Then
    RSize = UBound(MatrixMC, 1) - 1
ElseIf CurMatrixFlag = 11 Then
    RSize = UBound(MatrixL, 1) - 1
End If
Dim lTimerVal As Long
    If FirstPass = 1 Then TimerVal = 0
           If TimerVal > 0 Then
            lTimerVal = TimerVal - 40
           Else
            lTimerVal = 480
           End If
           DistD = RSize / MatZoom(CurMatrixFlag)
            ''PosS(0) = MatCoord(MatFlag, 0)
            'PosE(0) = PosS(0) + DistD
            'PosS(1) = MatCoord(MatFlag, 1)
            'PosE(1) = PosS(1) + DistD
            'XX = MatCoord(CurMatrixFlag, 0)
            XPicAddjust2 = (Form1.Picture26.ScaleHeight) / DistD
           Form1.Picture26.DrawMode = 7
           Dim Portion As Long, lPortion As Long
           Portion = 255 - 255 * (TimerVal ^ 0.5 / 520 ^ 0.5)
           lPortion = 255 - 255 * (lTimerVal ^ 0.5 / 520 ^ 0.5)
           'XPicAddjust2 = Picture26.ScaleHeight / UBound(MatrixMC, 1)
           Form1.Picture26.AutoRedraw = True
           'If CurMatrixFlag = 8 Then
               xcoord = (BPCoord(0) - MatCoord(CurMatrixFlag, 0)) * XPicAddjust2
               ycoord = (BPCoord(1) - MatCoord(CurMatrixFlag, 1)) * XPicAddjust2
            'ElseIf CurMatrixFlag = 11 Then
             '   xcoord = ((BPCoord(0) / Len(StrainSeq(0)) * RSize) - MatCoord(CurMatrixFlag, 0)) * XPicAddjust2
             '   ycoord = ((BPCoord(1) / Len(StrainSeq(0)) * RSize) - MatCoord(CurMatrixFlag, 1)) * XPicAddjust2
            'End If
           'Picture26.Refresh
           
           
           Form1.Picture26.DrawWidth = (Int(lTimerVal / 100) + 1)
           If FirstPass = 0 Then
            
                Form1.Picture26.Circle (xcoord, ycoord), lTimerVal / 40, RGB(lPortion, lPortion, lPortion)
                xcoord = (BPCoord(1) - MatCoord(CurMatrixFlag, 0)) * XPicAddjust2
                ycoord = (BPCoord(0) - MatCoord(CurMatrixFlag, 1)) * XPicAddjust2
                Form1.Picture26.Circle (xcoord, ycoord), lTimerVal / 40, RGB(lPortion, lPortion, lPortion)
            Else
                FirstPass = 0
                
           End If
           Form1.Picture26.DrawWidth = (Int(TimerVal / 100) + 1)
           xcoord = (BPCoord(0) - MatCoord(CurMatrixFlag, 0)) * XPicAddjust2
           ycoord = (BPCoord(1) - MatCoord(CurMatrixFlag, 1)) * XPicAddjust2
           Form1.Picture26.Circle (xcoord, ycoord), TimerVal / 40, RGB(Portion, Portion, Portion)
           xcoord = (BPCoord(1) - MatCoord(CurMatrixFlag, 0)) * XPicAddjust2
           ycoord = (BPCoord(0) - MatCoord(CurMatrixFlag, 1)) * XPicAddjust2
           Form1.Picture26.Circle (xcoord, ycoord), TimerVal / 40, RGB(Portion, Portion, Portion)
           'Form1.Picture26.AutoRedraw = True
           If Portion < 200 Then
                X = X
            End If
           Form1.Picture26.DrawWidth = 1
            Form1.Picture26.DrawMode = 13
End Sub
Public Sub TXover3()
    SS = GetTickCount
    AbortFlag = 0
    Dim TotMat2() As Double, PosCount As Long, CurrentPos As Long, X As Long, Y As Long
    Dim z As Long, B As Long, A As Long, OWinLen As Long, OStepSize As Long, MoveOver As Long
    Dim AF As Double, CF As Double, GF As Double, TF As Double, SumFirstDiff As Double, MeanFirstDiff As Double
    Dim ODir As String, BootName As String, Header As String, GetStringA As String, GetStringB As String
    Dim Boots() As String
    Dim SumVal As Double
    Dim DSS() As Double
    Dim FirstDiff() As Double
    Dim PlotPos As Long
    Dim Pict As Long
    Dim XFactor As Double, YScaleFactor As Double, PicHeight As Double
    Dim PntAPI As POINTAPI
    Dim LSeqs As Long
    Dim LTree As Long, StartT As Long, NumWins As Long, FF As Integer
    Dim RndNum2 As Long, RndNum As Long
    Dim NameLen As Integer, ttSeqNum() As Integer, Px() As Integer, xx1() As Integer, xx2() As Integer
    Dim WeightMod() As Long, Num1() As Long, Num2() As Long, Weight() As Long, Location() As Long, Ally() As Long, Alias() As Long
    Dim TMat2() As Double, Num() As Double, DistVal() As Double, Prod1() As Double, Prod2() As Double, Prod3() As Double, DEN() As Long
    Dim SHolder() As Byte
    Dim TreeOut As String
    Dim NodeOrder() As Long
    Dim NodeLen() As Double
    Dim DoneNode() As Long
    Dim TempNodeOrder() As Long
    Dim DstOut() As Integer, tSeqNum() As Integer
    Dim RootNode() As Long
    Dim PTOMat() As Double
'TOPower = 0
    LSeqs = Len(StrainSeq(0))
    Screen.MousePointer = 11
    Form1.Picture8.Enabled = False
    'Form1.SSPanel8.Enabled = False
    Form1.Picture10.Enabled = False
    Form1.Combo1.Enabled = False
    Form1.Command29.Enabled = False
    Form1.SSPanel2.Enabled = False
    Form1.Picture5.Enabled = False
    Form1.Command6.Enabled = False
    Form1.Command25.Enabled = True
    Form1.Picture7.Enabled = False
    TOPFlag = 0

    Dim TOSSFlag As Integer

    TOSSFlag = 0

    Dim TOXOSeq() As String

    ReDim TOXOSeq(Nextno)
    
    
    Dim TOXDiffPos() As Long, TOXPosDiff() As Long
    
    
    'Set up copies of sequences
    ReDim TOXDiffPos(LSeqs + 200)
    ReDim TOXPosDiff(LSeqs + 200)

    If TManFlag = 8 Then 'if it is a manual scan
        ReDim RevSeq(Nextno + 1)
        ReDim TOSeq(Nextno)
        ToNumSeqs = 0

        For X = 0 To Nextno

            If Selected(X) = 1 Then
                TOSeq(ToNumSeqs) = StrainSeq(X)
                RevSeq(ToNumSeqs) = X
                ToNumSeqs = ToNumSeqs + 1
            End If

        Next 'X

        ToNumSeqs = ToNumSeqs - 1
        ToNumSeqs = ToNumSeqs
        
    Else ' If it is a checking scan
        Dim TSeq(2) As String ', X As Long
        For X = 1 To Len(StrainSeq(0))
            'If X > 2050 Then
            '    X = X
            'End If
            If SeqNum(X, Seq1) <> 46 And SeqNum(X, Seq2) <> 46 And SeqNum(X, Seq3) <> 46 Then
                TSeq(0) = TSeq(0) + Chr(SeqNum(X, Seq1) - 1)
                TSeq(1) = TSeq(1) + Chr(SeqNum(X, Seq2) - 1)
                TSeq(2) = TSeq(2) + Chr(SeqNum(X, Seq3) - 1)
            Else
                TSeq(0) = TSeq(0) + "-"
                TSeq(1) = TSeq(1) + "-"
                TSeq(2) = TSeq(2) + "-"
            End If
        Next X
        If Seq1 = Seq2 Or Seq1 = Seq3 Or Seq2 = Seq3 Then
            ToNumSeqs = Nextno
            ReDim TOSeq(ToNumSeqs + 1)

            For X = 0 To ToNumSeqs
                TOSeq(X) = StrainSeq(X)
            Next 'X

        Else
            Form1.SSPanel1.Caption = "Constructing an Outlyer Sequence"
            ToNumSeqs = 3
            ReDim TOSeq(ToNumSeqs + 1)
            TOSeq(0) = TSeq(0)
            TOSeq(1) = TSeq(1)
            TOSeq(2) = TSeq(2)
            Count = 0
            Count2 = 0
            Y = 0

            For X = 1 To LSeqs
                TOXPosDiff(X) = Y

                If SeqNum(X, Seq1) = SeqNum(X, Seq2) And SeqNum(X, Seq1) = SeqNum(X, Seq3) Then
                    
                    TOSeq(3) = TOSeq(3) + Mid$(StrainSeq(TreeTraceSeqs(1, Seq1)), X, 1)
                Else
                    OddOne = -1
                    GoodOne = -1

                    If SeqNum(X, Seq1) = 46 Or SeqNum(X, Seq2) = 46 Or SeqNum(X, Seq3) = 46 Then
                        TOSeq(3) = TOSeq(3) + "-"
                        Mid$(TOSeq(0), X, 1) = "-"
                        Mid$(TOSeq(1), X, 1) = "-"
                        Mid$(TOSeq(2), X, 1) = "-"
                    ElseIf SeqNum(X, Seq1) <> SeqNum(X, Seq2) And SeqNum(X, Seq1) <> SeqNum(X, Seq3) And SeqNum(X, Seq2) <> SeqNum(X, Seq3) Then
                        TOSeq(3) = TOSeq(3) + "-"
                        Mid$(TOSeq(0), X, 1) = "-"
                        Mid$(TOSeq(1), X, 1) = "-"
                        Mid$(TOSeq(2), X, 1) = "-"
                    Else
                        Y = Y + 1
                        TOXDiffPos(Y) = X
                        TOXPosDiff(X) = Y
                        Y = Y + 1
                        TOXDiffPos(Y) = X
                        TOXPosDiff(X) = Y

                        If SeqNum(X, Seq1) = SeqNum(X, Seq2) Then
                            Mid$(TOSeq(0), X, 1) = "C"
                            Mid$(TOSeq(1), X, 1) = "C"
                            Mid$(TOSeq(2), X, 1) = "T"
                            TOSeq(3) = TOSeq(3) + "G"
                            TOXOSeq(0) = TOXOSeq(0) + "CA"
                            TOXOSeq(1) = TOXOSeq(1) + "CA"
                            TOXOSeq(2) = TOXOSeq(2) + "TA"
                            TOXOSeq(3) = TOXOSeq(3) + "GA"
                        ElseIf SeqNum(X, Seq1) = SeqNum(X, Seq3) Then
                            Mid$(TOSeq(0), X, 1) = "C"
                            Mid$(TOSeq(2), X, 1) = "C"
                            Mid$(TOSeq(1), X, 1) = "T"
                            TOSeq(3) = TOSeq(3) + "G"
                            TOXOSeq(0) = TOXOSeq(0) + "CA"
                            TOXOSeq(2) = TOXOSeq(2) + "CA"
                            TOXOSeq(1) = TOXOSeq(1) + "TA"
                            TOXOSeq(3) = TOXOSeq(3) + "GA"
                        ElseIf SeqNum(X, Seq2) = SeqNum(X, Seq3) Then
                            Mid$(TOSeq(2), X, 1) = "C"
                            Mid$(TOSeq(1), X, 1) = "C"
                            Mid$(TOSeq(0), X, 1) = "T"
                            TOSeq(3) = TOSeq(3) + "G"
                            TOXOSeq(2) = TOXOSeq(2) + "CA"
                            TOXOSeq(1) = TOXOSeq(1) + "CA"
                            TOXOSeq(0) = TOXOSeq(0) + "TA"
                            TOXOSeq(3) = TOXOSeq(3) + "GA"
                        End If

                    End If

                End If

                If X / 100 = CLng(X / 100) Then Form1.ProgressBar1.Value = (X / LSeqs) * 100
            Next 'X

        End If

    End If


    ReDim TOTSeqNum(Len(TOSeq(0)), ToNumSeqs)
    Dim TSeqSpaces() As Integer
    ReDim TSeqSpaces(Len(TOSeq(0)), ToNumSeqs)
    
    For X = 0 To ToNumSeqs
        Dummy = CopyString(Len(TOSeq(0)), TOTSeqNum(0, X), TOSeq(X), TSeqSpaces(0, X))
    Next X
    
    For X = 0 To Len(TOSeq(0))
        If TOTSeqNum(X, 0) = 46 Or TOTSeqNum(X, 1) = 46 Or TOTSeqNum(X, 2) = 46 Then
            TOTSeqNum(X, 0) = 46
            TOTSeqNum(X, 1) = 46
            TOTSeqNum(X, 2) = 46
            TOTSeqNum(X, 3) = 46
        End If
    Next X
    If TOSSFlag = 1 Then
        OWinLen = TOWinLen
        OStepSize = TOStepSize
        TOWinLen = 100
        TOStepSize = 10
        TOSeq(0) = TOXOSeq(0)
        TOSeq(1) = TOXOSeq(1)
        TOSeq(2) = TOXOSeq(2)
        TOSeq(3) = TOXOSeq(3)
    End If
    
    ReDim xx1(3)
    ReDim xx2(3)
    ReDim Prod1(Len(TOSeq(0)))
    ReDim Prod2(Len(TOSeq(0)))
    ReDim Prod3(Len(TOSeq(0)))
    ReDim Alias(Len(TOSeq(0)))
    ReDim Ally(Len(TOSeq(0)))
    ReDim Location(Len(TOSeq(0)))
    ReDim Px(ToNumSeqs, Len(TOSeq(0)))
    ReDim TTempSeq2(Len(TOSeq(0)), ToNumSeqs)
    ReDim DistVal(1)
    ReDim Num1(1)
    ReDim Num2(1)
    ReDim DEN(1)
    ReDim Num(1)
    ReDim WeightMod(0, Len(TOSeq(0)))
    ReDim TotMat(ToNumSeqs, ToNumSeqs)
    ReDim SHolder((ToNumSeqs) * 40 * 2)
    ReDim Weight(0, Len(TOSeq(0)))
    
    Form1.ProgressBar1.Value = 5
    Form1.SSPanel1.Caption = "Calculating Average Distance"
    
    ReDim LastMatrix(ToNumSeqs, ToNumSeqs)
    ReDim Boots(ToNumSeqs + 1)
    
    ODir = CurDir
    ChDir App.Path
    ChDrive App.Path
    'TOSeq(1) = X
    'Do Batch Files
    Open "dnadist.bat" For Output As #1
    Print #1, "dnadist <optfiled"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
Close #1
    
    Open "neighbor1.bat" For Output As #1
    Print #1, "neighborrdp <optfilen1"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
    Close #1
    
    Open "neighbor2.bat" For Output As #1
    Print #1, "neighborrdp <optfilen2"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
    Close #1
    Open "fitch1.bat" For Output As #1
    Print #1, "fitch <optfilef1"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
    Close #1
    Open "fitch2.bat" For Output As #1
    Print #1, "fitch <optfilef2"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
    Close #1
    Open "fitch3.bat" For Output As #1
    Print #1, "fitch <optfilef3"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
    Close #1
    Open "fitch4.bat" For Output As #1
    Print #1, "fitch <optfilef4"
    Print #1, "del outfilex"
    Print #1, "rename outfile outfilex"
    Close #1
    
    
    

    
    'Do Optionfiles
    
    'Neighbor option file - used during DSS calculations.
    Open "optfilen1" For Output As #1
    Print #1, "j"
    Print #1, TORndNum
    Print #1, "2"
    Print #1, "y"
    Close #1
    'Neighbor2 option file - used for parametric bootstrap
    Open "optfilen2" For Output As #1
    Print #1, "j"
    Print #1, TORndNum
    Print #1, "2"
    Print #1, "y"
    Close #1
    'Fitch1 option file
    Open "optfilef1" For Output As #1
    Print #1, "g"
    Print #1, "j"
    Print #1, TORndNum
    Print #1, "1"
    Print #1, "p"
    Print #1, TOPower
    Print #1, "2"
    Print #1, "y"
    Close #1
    'Fitch2 option file - used fo SS genration with "user"
    'defined tree
    Open "optfilef2" For Output As #1
    Print #1, "u"
    Print #1, "p"
    Print #1, TOPower
    Print #1, "2"
    Print #1, "y"
    Close #1
    'Fitch3 option file - used for parametric Bootstrap (LS)
    Open "optfilef3" For Output As #1
    Print #1, "g"
    Print #1, "j"
    Print #1, TORndNum
    Print #1, "3"
    Print #1, "p"
    Print #1, TOPower
    'Print #1, "2"
    Print #1, "y"
    Close #1
    'Fitch4 option file - used with NJ generated tree.
    Open "optfilef4" For Output As #1
    Print #1, "u"
    Print #1, "p"
    Print #1, TOPower
    Print #1, "2"
    Print #1, "y"
    Close #1
    'Clean up files
    On Error Resume Next
    Kill "outfile"
    Kill "outfilex"
    Kill "infile"
    Kill "dist.mat"
    Kill "distmatrix"
    Kill "simtree"
    Kill "infilebak"
    On Error GoTo 0
    Sleep (5)
    'STEP1 Work out average distances

    'MakeSubAlign NumberOfSeqs, Len(StrainSeq(0)), B, DPStep, DPWindow, TTempSeq2(0, 0), TempSeq2(0, 0)
    DNADIST TOCoeffVar, TOTvTs, TOFreqFlag, TOModel, TOFreqA, TOFreqC, TOFreqG, TOFreqT, ToNumSeqs + 1, Len(TOSeq(0)), TOTSeqNum(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), TotMat(0, 0)
    
    ReDim TotMat2(ToNumSeqs, ToNumSeqs)
            For z = 0 To ToNumSeqs - 1
                For Y = z + 1 To ToNumSeqs
                    'If TotMat(Z, Y) = 0 Then TotMat(Z, Y) = 0.00005
                    TotMat2(z, Y) = TotMat(z, Y)
                    TotMat2(Y, z) = TotMat(z, Y)
                Next Y
            Next z
    
    DoEvents
    
    If AbortFlag = 1 Then
        OldAbort = 1
        Form1.Command25.Enabled = False
        Form1.ProgressBar1.Value = 0
        Form1.SSPanel1.Caption = ""
        Form1.ProgressBar1.Value = 0
        Form1.SSPanel1.Caption = ""
        Form1.Picture10.Enabled = True
        Form1.Frame7.Enabled = True
        Form1.SSPanel2.Enabled = True
        Form1.Picture8.Enabled = True
        Form1.Picture5.Enabled = True
        Form1.Command6.Enabled = True
        On Error Resume Next
        Close #1
        Close #2
        On Error GoTo 0
        Exit Sub
    End If

    'Exit Sub

    'Call ReadDistMatrix(TONumSeqs, "outfilex")
    ReDim PTOMat(ToNumSeqs, ToNumSeqs)
    For A = 0 To ToNumSeqs
        For B = 0 To ToNumSeqs
            LastMatrix(A, B) = TotMat(A, B)
            PTOMat(A, B) = TotMat(A, B)
        Next 'B
    Next 'A

    Call AverageMatrix(ToNumSeqs, MatAverage)

    TAv = MatAverage
    'Name "outfilex" As "distmatrix"
    Sleep (5)
    ReDim SSScore(3, LSeqs)
    '2.5 seconds
    Form1.ProgressBar1 = 18
    'Get forward and backwards SSScores.

    Call TopalRun3(MatAverage, CurrentPos)

    If AbortFlag = 1 Then
        OldAbort = 1
        Form1.Command25.Enabled = False
        'ChDir (olddir$)
        Form1.ProgressBar1.Value = 0
        Form1.SSPanel1.Caption = ""
        Form1.ProgressBar1.Value = 0
        Form1.SSPanel1.Caption = ""
        'Form1.SSPanel8.Enabled = True
        Form1.Picture10.Enabled = True
        Form1.Frame7.Enabled = True
        Form1.SSPanel2.Enabled = True
        Form1.Picture8.Enabled = True
        Form1.Picture5.Enabled = True
        Form1.Command6.Enabled = True
        On Error Resume Next
        Close #1
        Close #2
        On Error GoTo 0
        Exit Sub
    End If
    Dim TPo As Long
    Dim TSm As Long

    TSm = TOSmooth - 1
    
    '14 seconds
    'Calculate DSS and Firstdiff
    Form1.ProgressBar1.Value = 100
    ReDim FirstDiff(CurrentPos)
    ReDim DSS(TOPerms, CurrentPos * 3)

    For X = 0 To CurrentPos - 1

        If Abs(SSScore(0, X) - SSScore(1, X)) > Abs(SSScore(2, X) - SSScore(3, X)) Then
            DSS(0, X) = Abs(SSScore(0, X) - SSScore(1, X))
        Else
            DSS(0, X) = Abs(SSScore(2, X) - SSScore(3, X))
        End If

        If X > 0 Then
            FirstDiff(X - 1) = DSS(0, X - 1) - DSS(0, X)
            SumFirstDiff = SumFirstDiff + FirstDiff(X - 1)
        End If

    Next 'X

    MeanFirstDiff = SumFirstDiff / (CurrentPos - 1)
    SumFirstDiff = 0

    For X = 0 To PosCount - 1
        SumFirstDiff = SumFirstDiff + (FirstDiff(X) - MeanFirstDiff) * (FirstDiff(X) - MeanFirstDiff)
    Next 'X

    VarFirstDiff = SumFirstDiff / (PosCount - 2)

    For X = 0 To CurrentPos - 1
        SumFirstDiff = SumFirstDiff + (FirstDiff(X) - MeanFirstDiff) * (FirstDiff(X) - MeanFirstDiff)
    Next 'X

    VarFirstDiff = SumFirstDiff / (CurrentPos - 2)
    'Look at "choosebig" to see the calculation of lower and upperbound
    '95% and 99% CIsI may need to write these to a file
    ' 95%CI = 1.96*(sqrt(varfirstdiff))
    ' 95%CI = 2.58*(sqrt(varfirstdiff))
    'Parametric bootstrapping bit
    '4.8 seconds
    TOEndPlot = CurrentPos / 2

    If TOPerms > 0 Then 'make simtree
        Form1.ProgressBar1.Value = 2
        Form1.SSPanel1.Caption = "Drawing LS Tree"
        'Draw the plot to show that something is happening

        If Seq1 = Seq2 Or Seq1 = Seq3 Or Seq2 = Seq3 Then
        Else
            Form1.Picture7.Top = 0
            Form1.Picture7.ScaleMode = 3
            Form1.Picture10.BackColor = BackColours
            Form1.Picture7.BackColor = BackColours
            Form1.Picture10.ScaleMode = 3
            YScaleFactor = 0.85
            PicHeight = Form1.Picture7.Height * YScaleFactor
            XFactor = ((Form1.Picture7.Width - 40) / LSeqs)
            
            
            SumVal = 0

            For z = 0 To TSm
                SumVal = SumVal + DSS(0, z)
            Next 'Z
            ReDim SmoothDSS(0, TOEndPlot * 2)
            SmoothDSS(0, 0) = SumVal / TOSmooth
            'Exit Sub
            For X = 1 To TOEndPlot 'CurrentPos - (TOSmooth - 1)
                SumVal = SumVal - DSS(0, X - 1) + DSS(0, X + TOSmooth - 1)
                SmoothDSS(0, X) = SumVal / TOSmooth
            Next 'X
            
            
            ReDim HighestDSS(TOPerms)

            For Y = 1 To TOEndPlot
                'For X = 0 To TOPerms

                If SmoothDSS(0, Y) >= HighestDSS(0) Then
                    HighestDSS(0) = SmoothDSS(0, Y)
                    
                End If

                'Next 'X
            Next 'Y

            TOHigh = HighestDSS(0)
            Form1.Picture10.BackColor = BackColours
            Form1.Picture7.BackColor = BackColours
            Form1.Picture10.ScaleMode = 3
            Form1.Picture7.ScaleMode = 3
            If ManFlag = -1 And TManFlag = -1 Then
            
                Call DoAxes(0, Len(StrainSeq(0)), -1, TOHigh, 0, 1, "DSS")
            Else
                Call DoAxes(0, Len(StrainSeq(0)), -1, TOHigh, 0, 0, "DSS")
            End If

            Form1.Picture7.Enabled = False
            PlotPos = Int((TOWinLen + TOStepSize * TSm) / 2)
            Form1.Picture7.DrawWidth = 3
            Form1.Picture7.ForeColor = RGB(180, 180, 180)
            Pict = Form1.Picture7.hDC
            
            If TOSSFlag = 0 Then
                Dummy = MoveToEx(Pict, 30 + (PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / HighestDSS(0)) * (PicHeight - 35)), PntAPI)
                'Exit Sub
                PlotPos = PlotPos + TOStepSize

                For X = 1 To TOEndPlot
                    Dummy = LineTo(Pict, 30 + PlotPos * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / HighestDSS(0)) * (PicHeight - 35)))
                    PlotPos = PlotPos + TOStepSize
                Next 'X

            Else
                Dummy = MoveToEx(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / HighestDSS(0)) * (PicHeight - 35)), PntAPI)
                PlotPos = PlotPos + TOStepSize
                X = 1

                For X = 1 To TOEndPlot
                    Dummy = LineTo(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / HighestDSS(0)) * (PicHeight - 35)))
                    PlotPos = PlotPos + TOStepSize
                Next 'X

            End If

            Dim EXCurrentPos As Long

            EXCurrentPos = CurrentPos
            Form1.Picture7.DrawWidth = 1
            Form1.Picture7.ForeColor = RGB(0, 0, 0)
            'Exit Sub
            Pict = Form1.Picture7.hDC
            PlotPos = Int((TOWinLen + TOStepSize * TSm) / 2)
            Form1.Picture7.ForeColor = RGB(0, 0, 0)

            If TOSSFlag = 0 Then
                Dummy = MoveToEx(Pict, 30 + (PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / HighestDSS(0)) * (PicHeight - 35)), PntAPI)
                PlotPos = PlotPos + TOStepSize

                For X = 1 To TOEndPlot
                    Dummy = LineTo(Pict, 30 + PlotPos * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / HighestDSS(0)) * (PicHeight - 35)))
                    PlotPos = PlotPos + TOStepSize
                Next 'X

            Else
                Dummy = MoveToEx(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / HighestDSS(0)) * (PicHeight - 35)), PntAPI)
                PlotPos = PlotPos + TOStepSize

                For X = 1 To TOEndPlot
                    Dummy = LineTo(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / HighestDSS(0)) * (PicHeight - 35)))
                    PlotPos = PlotPos + TOStepSize
                Next 'X

            End If

            If ManFlag = -1 And RunFlag >= 1 Then

                Call Highlight

            End If

            Form1.Picture7.Refresh
        End If

        On Error Resume Next
        Kill "infile"
        On Error GoTo 0
        
        Open "infile" For Output As #1
    Pos = 0
            Dim TS As String
            Print #1, "  " & CStr(ToNumSeqs + 1)

            For X = 0 To ToNumSeqs
                OutString = "S0" & Trim$(CStr(X)) & String$(9 - Len(Trim$(CStr(X))), " ")

                For Y = 0 To ToNumSeqs
                    
                    If PTOMat(X, Y) = 0 Then
                        OutString = OutString & "  0.00000"
                    ElseIf PTOMat(X, Y) = 1 Then
                        OutString = OutString & "  1.00000"
                    ElseIf PTOMat(X, Y) >= 10 Then
                        TS = left(Trim$(CStr(PTOMat(X, Y))), 7)
                        OutString = OutString & " " & TS
                        If Len(TS) < 7 Then
                            OutString = OutString & String$((7 - Len(TS)), "0")
                        End If
                    ElseIf PTOMat(X, Y) > 1 Then
                        TS = left(Trim$(CStr(PTOMat(X, Y))), 7)
                        OutString = OutString & "  " & TS
                        If Len(TS) < 7 Then
                            OutString = OutString & String$((7 - Len(TS)), "0")
                        End If
                    Else
                        TS = PTOMat(X, Y) * 10000
                        
                        
                        TS = CLng(PTOMat(X, Y) * 100000)
                        If Len(TS) <= 5 Then
                            TS = "0." & String(5 - Len(TS), "0") & TS
                        End If
                        'TS = "0." & String(5 - Len(TS), "0") & TS
                        TS = left(Trim$(TS), 7)
                        OutString = OutString & "  " & TS
                        If Len(TS) < 7 Then
                            OutString = OutString & String$((7 - Len(TS)), "0")
                        End If
                    End If
                Next 'Y
                Pos = InStr(1, OutString, ",", vbBinaryCompare)
                Do While Pos > 0
                    Mid(OutString, Pos, 1) = "."
                    Pos = InStr(1, OutString, ",", vbBinaryCompare)
                Loop
               
                Print #1, OutString
            Next 'x

        Close #1
        
        If TOTreeType = 1 Or ToNumSeqs = 2 Then
            
           
            Pos = 0
            
            'This is where the writing to file thing used to go
            ReDim SHolder((ToNumSeqs + 1) * 40 * 2)
            'Do NJ tree
            ReDim ColTotals(ToNumSeqs)
            
            ReDim TreeArray(ToNumSeqs, ToNumSeqs)
            LTree = NEIGHBOUR(1, 2, TORndNum, 1, ToNumSeqs + 1, TotMat2(0, 0), SHolder(0), ColTotals(0), TreeArray(0, 0))
            LTree = LTree + 2
            SHolder(LTree - 1) = 10
            SHolder(LTree) = 13
            f2 = FreeFile
            On Error Resume Next
            Kill "intree"
            On Error GoTo 0
            Open "intree" For Binary Access Write As #f2
            Put #f2, 1, SHolder
            Close #f2
            'ShellAndClose "neighbor2.bat", 0

            If AbortFlag = 1 Then
                OldAbort = 1
                Form1.Command25.Enabled = False
                'ChDir (olddir$)
                Form1.ProgressBar1.Value = 0
                Form1.SSPanel1.Caption = ""
                Form1.ProgressBar1.Value = 0
                Form1.SSPanel1.Caption = ""
                'Form1.SSPanel8.Enabled = True
                Form1.Picture10.Enabled = True
                Form1.Frame7.Enabled = True
                Form1.SSPanel2.Enabled = True
                Form1.Picture8.Enabled = True
                Form1.Picture5.Enabled = True
                Form1.Command6.Enabled = True
                On Error Resume Next
                Close #1
                Close #2
                On Error GoTo 0
                Exit Sub
            End If
            
            
        
    
            ShellAndClose "fitch4.bat", 0

            If AbortFlag = 1 Then
                OldAbort = 1
                Form1.Command25.Enabled = False
                'ChDir (olddir$)
                Form1.ProgressBar1.Value = 0
                Form1.SSPanel1.Caption = ""
                Form1.ProgressBar1.Value = 0
                Form1.SSPanel1.Caption = ""
                'Form1.SSPanel8.Enabled = True
                Form1.Picture10.Enabled = True
                Form1.Frame7.Enabled = True
                Form1.SSPanel2.Enabled = True
                Form1.Picture8.Enabled = True
                Form1.Picture5.Enabled = True
                Form1.Command6.Enabled = True
                On Error Resume Next
                Close #1
                Close #2
                On Error GoTo 0
                Exit Sub
            End If
            Name "outtree" As "simtree"
        Else
            On Error Resume Next
            Kill "outfile"
            Kill "outtree"
            On Error GoTo 0
            ShellAndClose "fitch3.bat", 0

            If AbortFlag = 1 Then
                OldAbort = 1
                Form1.Command25.Enabled = False
                'ChDir (olddir$)
                Form1.ProgressBar1.Value = 0
                Form1.SSPanel1.Caption = ""
                Form1.ProgressBar1.Value = 0
                Form1.SSPanel1.Caption = ""
                'Form1.SSPanel8.Enabled = True
                Form1.Picture10.Enabled = True
                Form1.Frame7.Enabled = True
                Form1.SSPanel2.Enabled = True
                Form1.Picture8.Enabled = True
                Form1.Picture5.Enabled = True
                Form1.Command6.Enabled = True
                On Error Resume Next
                Close #1
                Close #2
                On Error GoTo 0
                Exit Sub
            End If
            Name "outtree" As "simtree"
        End If

        On Error Resume Next
        Kill "distmatrix"
        On Error GoTo 0
        Name "infile" As "distmatrix"
        
        Sleep (5)
        TOPFlag = 1
        SS = GetTickCount
        Call TopalRunPerms3(MatAverage, CurrentPos)
        EE = GetTickCount
        TT = EE - SS '37
        '68 seconds

        If AbortFlag = 1 Then
            OldAbort = 1
            Form1.Command25.Enabled = False
            'ChDir (olddir$)
            Form1.ProgressBar1.Value = 0
            Form1.SSPanel1.Caption = ""
            Form1.ProgressBar1.Value = 0
            Form1.SSPanel1.Caption = ""
            'Form1.SSPanel8.Enabled = True
            Form1.Picture10.Enabled = True
            Form1.Frame7.Enabled = True
            Form1.SSPanel2.Enabled = True
            Form1.Picture8.Enabled = True
            Form1.Picture5.Enabled = True
            Form1.Command6.Enabled = True
            On Error Resume Next
            Close #1
            Close #2
            On Error GoTo 0
            Exit Sub
        End If

        Form1.ProgressBar1.Value = 100
        X = 1
        Y = 0
        ttCPos = CurrentPos
        CurrentPos = TOEndPlot 'CLng((CurrentPos / 2)) ' / TOPerms)

        Do

            If Abs(SSScore(0, MoveOver + Y) - SSScore(1, MoveOver + Y)) > Abs(SSScore(2, MoveOver + Y) - SSScore(3, MoveOver + Y)) Then
                DSS(X, Y) = Abs(SSScore(0, MoveOver + Y) - SSScore(1, MoveOver + Y))
            Else
                DSS(X, Y) = Abs(SSScore(2, MoveOver + Y) - SSScore(3, MoveOver + Y))
                'GoTo FinishAndCleanUp
            End If

            Y = Y + 1

            If Y >= CurrentPos Or Y > UBound(DSS, 2) Or Y + MoveOver > UBound(SSScore, 2) Then
                MoveOver = MoveOver + Y
                Y = 0
                X = X + 1
                'MoveOver = (X - 1) * CurrentPos

                If X > TOPerms Or Y + MoveOver > UBound(SSScore, 2) Then Exit Do
            End If

        Loop

    Else
        CurrentPos = CLng(CurrentPos / 2)
    End If

    If EXCurrentPos > 0 Then EXCurrentPos = CLng(EXCurrentPos / 2)

    If EXCurrentPos < CurrentPos Then
        ReDim SmoothDSS(TOPerms, TOEndPlot * 2)
    Else
        ReDim SmoothDSS(TOPerms, TOEndPlot * 2)
    End If

    ReDim HighestDSS(TOPerms)
    'Smooth DSS values

    Dim TCurrentPos As Long
    

    TPo = TOEndPlot - (TOSmooth - 1)

    
    'If TOSmooth > 1 Then

    For Y = 0 To TOPerms

        If Y = 0 And EXCurrentPos > 0 Then
            TCurrentPos = CurrentPos
            CurrentPos = EXCurrentPos
        End If

        SumVal = 0

        For z = 0 To TSm
            SumVal = SumVal + DSS(Y, z)
        Next 'Z

        SmoothDSS(Y, 0) = SumVal / TOSmooth

        For X = 1 To TOEndPlot 'CurrentPos - (TOSmooth - 1)
            SumVal = SumVal - DSS(Y, X - 1) + DSS(Y, X + TOSmooth - 1)
            SmoothDSS(Y, X) = SumVal / TOSmooth
        Next 'X

        If TCurrentPos > 0 Then CurrentPos = TCurrentPos
    Next 'Y

    'End If
    'For X = 0 To CurrentPos
    '    If SmoothDSS(0, X) <> DSS(0, X) Then
    '
    '    End If
    'Next 'X
    ReDim HighestDSS(TOPerms)

    For X = 0 To TOPerms

        If X = 0 And EXCurrentPos > 0 Then
            TCurrentPos = CurrentPos
            CurrentPos = EXCurrentPos
        End If

        For Y = 0 To TPo 'CurrentPos - (TOSmooth - 1)

            If SmoothDSS(X, Y) >= HighestDSS(X) Then
                HighestDSS(X) = SmoothDSS(X, Y)
            End If

        Next 'Y

        If TCurrentPos > 0 Then CurrentPos = TCurrentPos
    Next 'X

    TOHigh = HighestDSS(0)

    For X = 1 To TOPerms

        If HighestDSS(X) > TOHigh Then TOHigh = HighestDSS(X)
    Next 'X

    '(2)for every permutation take dss values and calculate
    'moving average along the sequence (with a set window length)
    'write highest average to an array(x).
    '(3) sort the set of highest values from lowest to highest.
    '(4)from actual dataset .dss values smooth the values (moving
    'average etc etc.)
    '(5) significance = proportion of simulated .dss values that are
    'higher than the calculated values. eg-pvalue cutoff =0.5 then the
    '0.05 dss cutoff = the 5th percentile of simulated .dss values
    ChDir ODir
    ChDrive ODir

    If TOHigh = 0 Then
        GoTo FinishAndCleanUp:
    End If
    
    Dim MaxDSS() As Double
    ReDim MaxDSS(TOPerms)
    For X = 1 To TOPerms
        For Y = 1 To TPo
            If MaxDSS(X) < SmoothDSS(X, Y) Then MaxDSS(X) = SmoothDSS(X, Y)
        Next Y
    Next X
    'sort the values
    For X = 1 To TOPerms
        Win = -1: winp = -1
        For Y = X To TOPerms
            If MaxDSS(Y) > Win Then
                Win = MaxDSS(Y)
                winp = Y
            End If
        Next Y
        MaxDSS(winp) = MaxDSS(X)
        MaxDSS(X) = Win
        
    Next X
    ReDim GCritval(10)
    GCritval(2) = MaxDSS(CLng(0.01 * TOPerms + 0.5))
    GCritval(3) = MaxDSS(CLng(0.05 * TOPerms + 0.5))
    
    Form1.ProgressBar1.Value = 100

    If (Seq1 = Seq2 Or Seq1 = Seq3 Or Seq2 = Seq3) And ManFlag = -1 And TManFlag <> 8 Then
    Else
        Form1.Picture7.Top = 0
        Form1.Picture7.ScaleMode = 3
        Form1.Picture10.BackColor = BackColours
        Form1.Picture7.BackColor = BackColours
        Form1.Picture10.ScaleMode = 3
        XFactor = ((Form1.Picture7.Width - 40) / LSeqs)
        YScaleFactor = 0.85
        PicHeight = Form1.Picture7.Height * YScaleFactor

        If ManFlag = -1 And TManFlag = -1 Then
                Call FindSubSeqBS
                Call DoAxes(0, Len(StrainSeq(0)), -1, TOHigh, 0, 1, "DSS")
            Else
                
                Call DoAxes(0, Len(StrainSeq(0)), -1, TOHigh, 0, 0, "DSS")
            End If

        If EXCurrentPos > 0 Then CurrentPos = EXCurrentPos
        Form1.Picture7.DrawWidth = 3
        Form1.Picture7.ForeColor = RGB(180, 180, 180)
        Pict = Form1.Picture7.hDC
        PlotPos = Int((TOWinLen + TOStepSize * TSm) / 2)
        If TOSSFlag = 0 Then
            Dummy = MoveToEx(Pict, (30 + PlotPos * XFactor), PicHeight - (15 + (SmoothDSS(0, 0) / TOHigh) * (PicHeight - 35)), PntAPI)
            'Exit Sub
            PlotPos = PlotPos + TOStepSize
            Form1.Picture7.AutoRedraw = True

            For X = 1 To TPo
                Dummy = LineTo(Pict, 30 + PlotPos * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / TOHigh) * (PicHeight - 35)))
                PlotPos = PlotPos + TOStepSize
            Next 'X

        Else
            Dummy = MoveToEx(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / TOHigh) * (PicHeight - 35)), PntAPI)
            PlotPos = PlotPos + TOStepSize

            For X = 1 To TPo
                Dummy = LineTo(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / TOHigh) * (PicHeight - 35)))
                PlotPos = PlotPos + TOStepSize
            Next 'X

        End If

        Form1.Picture7.DrawWidth = 1
        Form1.Picture7.ForeColor = RGB(0, 0, 0)

        If TCurrentPos > 0 And TCurrentPos >= CurrentPos * 2 Then
            CurrentPos = TCurrentPos / 2
        ElseIf TCurrentPos > 0 Then
            CurrentPos = TCurrentPos
        End If

        Form1.Picture7.ForeColor = QuaterColour

        Dim PFactor As Double

        If LSeqs < 32000 Then
            PFactor = XFactor
        Else
            PFactor = ((Form1.Picture7.Width - 40) / 32000)
        End If

        For Y = 1 To TOPerms
            PlotPos = Int((TOWinLen + TOStepSize * TSm) / 2)
            

            If TOSSFlag = 0 Then
                'Form1.Picture7.ForeColor = QBColor(4)
                Pict = Form1.Picture7.hDC
                Dummy = MoveToEx(Pict, 30 + PlotPos * PFactor, PicHeight - (15 + (SmoothDSS(Y, 0) / TOHigh) * (PicHeight - 35)), PntAPI)

                For X = 1 To TPo
                    PlotPos = PlotPos + TOStepSize
                    Dummy = LineTo(Pict, 30 + PlotPos * PFactor, PicHeight - (15 + (SmoothDSS(Y, X) / TOHigh) * (PicHeight - 35)))
                Next 'X

            Else
                Dummy = MoveToEx(Pict, 30 + TOXDiffPos(PlotPos) * PFactor, PicHeight - (15 + (SmoothDSS(Y, 0) / TOHigh) * (PicHeight - 35)), PntAPI)

                For X = 1 To TPo
                    PlotPos = PlotPos + TOStepSize
                    Dummy = LineTo(Pict, 30 + TOXDiffPos(PlotPos) * PFactor, PicHeight - (15 + (SmoothDSS(Y, X) / TOHigh) * (PicHeight - 35)))
                Next

            End If

        Next
        
        'critical vals
        Form1.Picture7.DrawStyle = 2
        Form1.Picture7.ForeColor = 0
        Dummy = MoveToEx(Pict, 30, PicHeight - (15 + (GCritval(3) / TOHigh) * (PicHeight - 35)), PntAPI)
        Dummy = LineTo(Pict, 30 + Len(StrainSeq(0)), PicHeight - (15 + (GCritval(3) / TOHigh) * (PicHeight - 35)))
        
        Dummy = MoveToEx(Pict, 30, PicHeight - (15 + (GCritval(2) / TOHigh) * (PicHeight - 35)), PntAPI)
        Dummy = LineTo(Pict, 30 + Len(StrainSeq(0)), PicHeight - (15 + (GCritval(2) / TOHigh) * (PicHeight - 35)))
        Form1.Picture7.DrawStyle = 0
        PlotPos = Int((TOWinLen + TOStepSize * TSm) / 2)

        If EXCurrentPos > 0 Then CurrentPos = TCurrentPos
        Form1.Picture7.ForeColor = RGB(0, 0, 0)

        If TOSSFlag = 0 Then
            Dummy = MoveToEx(Pict, 30 + PlotPos * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / TOHigh) * (PicHeight - 35)), PntAPI)
            PlotPos = PlotPos + TOStepSize

            For X = 1 To TPo
                Dummy = LineTo(Pict, 30 + PlotPos * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / TOHigh) * (PicHeight - 35)))
                PlotPos = PlotPos + TOStepSize
            Next

        Else
            Dummy = MoveToEx(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, 0) / TOHigh) * (PicHeight - 35)), PntAPI)
            PlotPos = PlotPos + TOStepSize

            For X = 1 To TPo
                Dummy = LineTo(Pict, 30 + TOXDiffPos(PlotPos) * XFactor, PicHeight - (15 + (SmoothDSS(0, X) / TOHigh) * (PicHeight - 35)))
                PlotPos = PlotPos + TOStepSize
            Next

        End If

        'Next 'X






        If ManFlag = -1 And RunFlag >= 1 Then

            Call Highlight

        End If

        Form1.Picture7.Refresh
    End If

        PosCount = TPo
'Get everything into the standard format for printing and saving
        NumLines = TOPerms + 1 'number of lines to print
        ReDim GPrint(NumLines - 1, PosCount + 1), GPrintCol(NumLines - 1), GPrintPos(NumLines - 1, PosCount + 1)
        If ManFlag = -1 And TManFlag = -1 Then
            ReDim GVarPos(0, LenXoverSeq)
            For X = 1 To LenXoverSeq
                GVarPos(0, X) = XDiffpos(X)
            Next X
        Else
            ReDim GVarPos(0, 0)
        End If
        
        
        GLegend = "DSS"
        GPrintLen = PosCount + 1 'how many points to plot
        GPrintCol(0) = 0 'line is yellow
        For X = 1 To NumLines - 1
            GPrintCol(X) = RGB(128, 128, 128)
        Next X
        GPrintNum = NumLines - 1 'six lines
        GPrintType = 0 'a normal line plot
        
        GPrintMin(0) = 0 'bottom val
        GPrintMin(1) = TOHigh 'upper val
        
        
        For Y = 0 To NumLines - 1
            PlotPos = Int((TOWinLen + TOStepSize * TSm) / 2)
            For X = 1 To PosCount
                
                GPrint(Y, X) = SmoothDSS(Y, X)
                GPrintPos(Y, X) = PlotPos
                PlotPos = PlotPos + TOStepSize
                
            Next X
        Next Y
        
        For X = 0 To NumLines - 1
            GPrintPos(X, GPrintLen) = Len(StrainSeq(0))
            GPrint(X, GPrintLen) = (GPrint(X, GPrintLen - 1) + GPrint(X, 1)) / 2
            GPrintPos(X, 0) = 1
            GPrint(X, 0) = (GPrint(X, GPrintLen - 1) + GPrint(X, 1)) / 2
        Next X
        GExtraTNum = 1
        ReDim GExtraText(GExtraTNum)
        GExtraText(0) = "Real sequences"
        GExtraText(1) = "Sequences simulated without recombination"
        If ManFlag = -1 Then
            GExtraText(0) = GExtraText(0) & " (" + StraiName(Seq1) + ", " + StraiName(Seq2) + ", " + StraiName(Seq3) + ", and a simulated outlyer)"
        Else
            GBlockNum = -1
        End If
        

EE = GetTickCount
TT = EE - SS



'46.6 -50 reps, kom, ns ,reu ,tas
'13.532
    'If TOEndPlot <= 2 Then TOEndPlot = PlotPos - TOStepSize
FinishAndCleanUp:
    Form1.ProgressBar1.Value = 0
    Form1.SSPanel1.Caption = ""
    Form1.Combo1.Enabled = True
    'Form1.SSPanel8.Enabled = True
    Form1.Picture10.Enabled = True
    Form1.Frame7.Enabled = True
    Form1.SSPanel2.Enabled = True
    Form1.Picture8.Enabled = True
    Form1.Picture5.Enabled = True
    Form1.Command6.Enabled = True
    Form1.Command25.Enabled = False
    Form1.Picture7.Enabled = True
    Form1.Command29.Enabled = True
    
    
    Form1.Picture21.PaintPicture Form1.Picture7.Image, Form1.Picture7.left, Form1.Picture7.Top + 5

    If TOSSFlag = 1 Then
        TOWinLen = OWinLen
        TOStepSize = OStepSize
    End If

End Sub
Public Sub TopalRun3(MatAveragex As Double, CPos As Long)
Dim NextWrite As Long
Dim OutString As String
    Dim NumDatasets As Long, C As Long, d As Long
    Dim NewAverage As Double, MultFactorx As Double
    Dim StartPos As Long, X As Long, Y As Long, z As Long, A As Long, B As Long
    Dim BootName As String
    Dim GetString As String
    Dim GetStringA As String
    Dim GetStringB As String
    Dim GetStringC As String
    Dim GetStringD As String
    Dim GetStringE As String
    Dim GetStringF As String
    Dim Boots() As String
    Dim MatrixByte() As Byte
    Dim Pos As Long
    Dim MatrixLen As Long, SeqNumber As Long
    Dim Len1 As Long
    Dim LastPos As Long
    Dim ReadposA As Long
    Dim LenString As Long
     Dim LTree As Long, StartT As Long, NumWins As Long
    Dim RndNum2 As Long, RndNum As Long
    Dim NameLen As Integer, ttSeqNum() As Integer, Px() As Integer, xx1() As Integer, xx2() As Integer
    Dim WeightMod() As Long, Num1() As Long, Num2() As Long, Weight() As Long, Location() As Long, Ally() As Long, Alias() As Long
    Dim TMat2() As Double, Num() As Double, DistVal() As Double, Prod1() As Double, Prod2() As Double, Prod3() As Double, DEN() As Long
    Dim SHolder() As Byte
    Dim TreeOut As String
    Dim NodeOrder() As Long
    Dim NodeLen() As Double
    Dim DoneNode() As Long
    Dim TempNodeOrder() As Long
    Dim DstOut() As Integer, tSeqNum() As Integer
    Dim RootNode() As Long
    
    
    ReDim xx1(3)
    ReDim xx2(3)
    ReDim Prod1(TOWinLen)
    ReDim Prod2(TOWinLen)
    ReDim Prod3(TOWinLen)
    ReDim Alias(TOWinLen)
    ReDim Ally(TOWinLen)
    ReDim Location(TOWinLen)
    ReDim Px(ToNumSeqs, TOWinLen)
    ReDim TTempSeq2(TOWinLen, ToNumSeqs)
    ReDim DistVal(1)
    ReDim Num1(1)
    ReDim Num2(1)
    ReDim DEN(1)
    ReDim Num(1)
    ReDim WeightMod(0, Int(TOWinLen / 2))
    ReDim TotMat(ToNumSeqs, ToNumSeqs)
    ReDim SHolder((ToNumSeqs) * 40 * 2)
    ReDim Weight(0, Int(TOWinLen / 2))
     SeqLen = Len(StrainSeq(0))

    ReDim Boots(ToNumSeqs + 1)
    
    StartPos = 1
    CPos = 0
   
    Form1.ProgressBar1.Value = 20

    'Now can sort out the optfiles

    

    Form1.SSPanel1.Caption = "Calculating Distances"

    If AbortFlag = 1 Then
        Exit Sub
    End If

    Form1.ProgressBar1.Value = 40
    StartPos = 1
    CPos = 0
    StartPos = 0
     Dim TempSeq2() As Integer
    ReDim TempSeq2(Int(TOWinLen / 2), ToNumSeqs)
    ReDim TotMat(ToNumSeqs, ToNumSeqs)
    Dim TS As String
    
  If TOTreeType = 1 Then
    f2 = FreeFile
        On Error Resume Next
        Kill "intree"
        On Error GoTo 0
        Open "intree" For Binary Access Write As #f2
        
    End If
    FF = FreeFile
    NumDatasets = 0
    NextWrite = 1
    Open "infile" For Output As #FF
        Do While StartPos + TOWinLen < SeqLen
            
            'Do the first half of the window.
            
            MakeSubAlign ToNumSeqs, Len(TOSeq(0)), StartPos + Int(TOWinLen / 4), 1, Int(TOWinLen / 2), TempSeq2(0, 0), TOTSeqNum(0, 0)
            
            DNADIST TOCoeffVar, TOTvTs, TOFreqFlag, TOModel, TOFreqA, TOFreqC, TOFreqG, TOFreqT, ToNumSeqs + 1, Int(TOWinLen / 2), TempSeq2(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), TotMat(0, 0)
            
            Call AverageMatrix(ToNumSeqs, NewAverage)
            
            If NewAverage > 0 Then
                MultFactorx = MatAveragex / NewAverage

                Call MatrixMultiply(ToNumSeqs, MultFactorx)

                For C = 0 To ToNumSeqs

                    For d = 0 To ToNumSeqs
                        LastMatrix(C, d) = TotMat(C, d)
                    Next 'd

                Next 'C

            Else
                MultFactorx = 1

                For C = 0 To ToNumSeqs

                    For d = 0 To ToNumSeqs
                        TotMat(C, d) = LastMatrix(C, d)
                    Next 'd

                Next 'C

            End If
            
            Call DistsToFile(FF, ToNumSeqs, TotMat())
            
            If TOTreeType = 1 Then
                ReDim TotMat2(ToNumSeqs, ToNumSeqs)
                For z = 0 To ToNumSeqs - 1
                    For Y = z + 1 To ToNumSeqs
                        'If TotMat(Z, Y) = 0 Then TotMat(Z, Y) = 0.00005
                        TotMat2(z, Y) = TotMat(z, Y)
                        TotMat2(Y, z) = TotMat(z, Y)
                    Next Y
                Next z
                
                Pos = 0
                
                'This is where the writing to file thing used to go
                ReDim SHolder((ToNumSeqs + 1) * 40 * 2)
                'Do NJ tree
                ReDim ColTotals(ToNumSeqs)
                
                ReDim TreeArray(ToNumSeqs, ToNumSeqs)
                LTree = NEIGHBOUR(1, 2, TORndNum, 1, ToNumSeqs + 1, TotMat(0, 0), SHolder(0), ColTotals(0), TreeArray(0, 0))
                LTree = LTree + 2
                SHolder(LTree - 1) = 10
                SHolder(LTree) = 13
                Put #f2, NextWrite, SHolder
                
                NextWrite = NextWrite + LTree
                
            End If
            
            'Do second half of window
            MakeSubAlign ToNumSeqs, Len(TOSeq(0)), StartPos + Int(TOWinLen / 2) + Int(TOWinLen / 4), 1, Int(TOWinLen / 2), TempSeq2(0, 0), TOTSeqNum(0, 0)
            For X = 0 To TOWinLen / 2
                TOTSeqNum(X + 100, 0) = TOTSeqNum(X + 100, 0)
                TempSeq2(X, 0) = TempSeq2(X, 0)
                TempSeq2(X, 1) = TempSeq2(X, 1)
                TempSeq2(X, 2) = TempSeq2(X, 2)
                TempSeq2(X, 3) = TempSeq2(X, 3)
            Next X
            DNADIST TOCoeffVar, TOTvTs, TOFreqFlag, TOModel, TOFreqA, TOFreqC, TOFreqG, TOFreqT, ToNumSeqs + 1, Int(TOWinLen / 2), TempSeq2(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), TotMat(0, 0)
            
            Call AverageMatrix(ToNumSeqs, NewAverage)

            If NewAverage > 0 Then
                MultFactorx = MatAveragex / NewAverage

                Call MatrixMultiply(ToNumSeqs, MultFactorx)

                For C = 0 To ToNumSeqs

                    For d = 0 To ToNumSeqs
                        LastMatrix(C, d) = TotMat(C, d)
                    Next 'd

                Next 'C

            Else
                MultFactorx = 1

                For C = 0 To ToNumSeqs

                    For d = 0 To ToNumSeqs
                        TotMat(C, d) = LastMatrix(C, d)
                    Next 'd

                Next 'C

            End If
            
            Call DistsToFile(FF, ToNumSeqs, TotMat())
            If TOTreeType = 1 Then
                ReDim TotMat2(ToNumSeqs, ToNumSeqs)
                For z = 0 To ToNumSeqs - 1
                    For Y = z + 1 To ToNumSeqs
                        'If TotMat(Z, Y) = 0 Then TotMat(Z, Y) = 0.00005
                        TotMat2(z, Y) = TotMat(z, Y)
                        TotMat2(Y, z) = TotMat(z, Y)
                    Next Y
                Next z
                
                Pos = 0
                
                'This is where the writing to file thing used to go
                ReDim SHolder((ToNumSeqs + 1) * 40 * 2)
                'Do NJ tree
                ReDim ColTotals(ToNumSeqs)
                
                ReDim TreeArray(ToNumSeqs, ToNumSeqs)
                LTree = NEIGHBOUR(1, 2, TORndNum, 1, ToNumSeqs + 1, TotMat(0, 0), SHolder(0), ColTotals(0), TreeArray(0, 0))
                LTree = LTree + 2
                SHolder(LTree - 1) = 10
                SHolder(LTree) = 13
                Put #f2, NextWrite, SHolder
                
                NextWrite = NextWrite + LTree
                
            End If
            CPos = CPos + 1
            StartPos = StartPos + TOStepSize
            NumDatasets = NumDatasets + 2
        Loop

    

    Close #FF
    If TOTreeType = 1 Then
        Close #f2
    End If
    StartPos = 1
    Form1.ProgressBar1.Value = 44
    '3.7 seconds
    CPos = 0
    Form1.SSPanel1.Caption = "Calculating SS Scores"
    BB = GetTickCount
    
    If TOPFlag = 0 Then
        
        'Fitch1 option file
        Open "optfilef1" For Output As #1
        Print #1, "g"
        Print #1, "j"
        Print #1, TORndNum
        Print #1, "1"
        Print #1, "p"
        Print #1, TOPower
        Print #1, "m"
        Print #1, NumDatasets
        Print #1, "3"
        Print #1, "2"
        Print #1, "3"
        Print #1, "y"
        Close #1
        'Fitch2 option file - used with NJ generated tree
        'and SS genration with "user" defined tree
        Open "optfilef2" For Output As #1
        Print #1, "u"
        Print #1, "p"
        Print #1, TOPower
        Print #1, "m"
        Print #1, NumDatasets
        Print #1, "3"
        Print #1, "2"
        Print #1, "3"
        Print #1, "y"
        Close #1
        Open "optfilen1" For Output As #1
        Print #1, "j"
        Print #1, TORndNum
        Print #1, "2"
        Print #1, "m"
        Print #1, NumDatasets
        'Print #1, "2"
        'Print #1, "3"
        Print #1, "3"
        Print #1, "y"
        Close #1
    End If
    
    If TOTreeType = 0 Then 'do NJ if it is either set or you're working with permutations.
       ' Do
        If X = 12345 Then
            ShellAndClose "neighbor1.bat", 0
            If AbortFlag = 1 Then
                AbortFlag = 1
                Exit Sub
            End If
            Open "treefile" For Binary Access Read As #1
    
            If LOF(1) = 0 Then
                MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component NeighborRDP.")
                AbortFlag = 1
                Close #1
                Exit Sub
            End If
        
            
            Close #1
            Kill "intree"
            Name "treefile" As "intree"
        End If
        FileCopy "infile", "infilebak"
        
        
        'Loop
        
        On Error Resume Next
            Kill "outfile"
            Kill "outtree"
        On Error GoTo 0
        
        
        
        
        Open "optfilef2" For Output As #1
        Print #1, "u"
        Print #1, "p"
        Print #1, TOPower
        Print #1, "m"
        Print #1, NumDatasets
        Print #1, "3"
        Print #1, "3"
        Print #1, "2"
        Print #1, "y"
        Close #1
        
        
        
        ShellAndClose "fitch2.bat", 0
        
        If AbortFlag = 1 Then
            Exit Sub
        End If
        On Error Resume Next
        Kill "infile"
        On Error GoTo 0
        Name "infilebak" As "infile"
    Else
        
        ShellAndClose "fitch1.bat", 0

        If AbortFlag = 1 Then
            Exit Sub
        End If

    End If

    Form1.ProgressBar1.Value = 78
    '4.9 seconds (NJ)
    '2.4 seconds (LS)
    Open "outfilex" For Binary Access Read As #1

    If LOF(1) = 0 Then
        MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component FITCH.")
        AbortFlag = 1
        Close #1
        Exit Sub
    End If
    
    GetStringA = String$(LOF(1), " ")
    Get #1, 1, GetStringA
    Close #1
    'get ssscores
    ReDim SSScore(3, NumDatasets)
    Pos = 1
    CPos = 0

    Do
        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            On Error Resume Next
            SSScore(0, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(0, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
        End If

        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            On Error Resume Next
            SSScore(2, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(2, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
            CPos = CPos + 1
        End If

    Loop

    Form1.SSPanel1.Caption = "Recalculating SS Scores"
    Open "infile" For Binary Access Read As #1

    If LOF(1) = 0 Then
        MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component FITCH.")
        AbortFlag = 1
        Close #1
        Exit Sub
    End If

    GetStringA = String$(LOF(1), " ")
    Get #1, 1, GetStringA
    Close #1
    
    'Work out distance matrix spacing
    Pos = InStr(1, GetStringA, "S00", vbBinaryCompare)
    Form1.Refresh
    HeaderLen = Pos - 1
    LastPos = Pos
    Pos = InStr(LastPos + 1, GetStringA, "S00", vbBinaryCompare)
    Len1 = Pos - LastPos
    
    LenString = Len(GetStringA)
    Pos = 1
    LastPos = -2
    ReadposA = 1
    Open "infile" For Output As #1
    
    NumDatasets = 0
    Do While (ReadposA + Len1 * 2) - 1 <= LenString
        GetStringC = Mid$(GetStringA, ReadposA, Len1)
        ReadposA = ReadposA + Len1
        GetStringE = Mid$(GetStringA, ReadposA, Len1)
        ReadposA = ReadposA + Len1
        Pos = InStr(LastPos + 3, GetStringB, ";", vbBinaryCompare)

        
        LastPos = Pos
        Print #1, GetStringE
        
        Print #1, GetStringC
        
        NumDatasets = NumDatasets + 2
    Loop
    
    Close #1

        
    ShellAndClose "fitch2.bat", 0
    
    If AbortFlag = 1 Then
        Exit Sub
    End If

    Form1.ProgressBar1.Value = 95

    Open "outfilex" For Binary Access Read As #1

    If LOF(1) = 0 Then
        MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component FITCH.")
        AbortFlag = 1
        Close #1
        Exit Sub
    End If

    GetStringA = String$(LOF(1), " ")
    Get #1, 1, GetStringA
    Close #1
    'get ssscores
    Pos = 1
    CPos = 0

    Do
        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            On Error Resume Next
            SSScore(1, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(1, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
        End If

        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            On Error Resume Next
            SSScore(3, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(3, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
            CPos = CPos + 1
        End If

    Loop

    CPos = NumDatasets
    'SSScore(3, 15) = X
End Sub
Public Sub DistsToFile(FF, ToNumSeqs, TotMat() As Double)

Dim Pos As Long, z As Long, Y As Long
Dim OutString As String, TS As String

Pos = 0
Print #FF, ""
Print #FF, "  " & CStr(ToNumSeqs + 1)
For z = 0 To ToNumSeqs
    If z < 10 Then
        OutString = "S0" & Trim$(CStr(z)) & String$(9 - Len(Trim$(CStr(z))), " ")
    Else
        OutString = "S" & Trim$(CStr(z)) & String$(9 - Len(Trim$(CStr(z))), " ")
    End If
    For Y = 0 To ToNumSeqs
        If TotMat(z, Y) = 0 Then
            OutString = OutString & "  0.0000"
        ElseIf TotMat(z, Y) = 1 Then
            OutString = OutString & "  1.0000"
        ElseIf TotMat(z, Y) >= 10 Then
            TS = left(Trim$(CStr(TotMat(z, Y))), 6)
            OutString = OutString & " " & TS
            If Len(TS) < 7 Then
                OutString = OutString & String$((7 - Len(TS)), "0")
            End If
        ElseIf TotMat(z, Y) > 1 Then
            TS = left(Trim$(CStr(TotMat(z, Y))), 5)
            OutString = OutString & "  " & TS
            If Len(TS) < 6 Then
                OutString = OutString & String$((6 - Len(TS)), "0")
            End If
        Else
            TS = left(Trim$(CStr(TotMat(z, Y))), 5)
            OutString = OutString & "  " & TS
            If Len(TS) < 6 Then
                OutString = OutString & String$((6 - Len(TS)), "0")
            End If
        End If
    Next 'Y
    Pos = InStr(1, OutString, ",", vbBinaryCompare)
    Do While Pos > 0
        Mid(OutString, Pos, 1) = "."
        Pos = InStr(1, OutString, ",", vbBinaryCompare)
    Loop
    Print #FF, OutString
Next 'z
End Sub

Public Sub TopalRunPerms3(MatAveragex As Double, CPos As Long)
    'This code will take a sequence (in TOSeq), split it into window
    'size chunks, calculate SS scores for chunks and return
    'an array containing the SS data
    
    Dim AF As Double, CF As Double, GF As Double, TF As Double, FF As Integer, FF2 As Integer
    Dim Target As String, OutString As String, GetString As String, GetStringA As String, GetStringB As String, GetStringC As String, GetStringD As String, GetStringE As String, GetStringF As String, BootName As String, Header As String
    Dim StartPos  As Long, SeqLen As Long, NumDatasets As Long, Pos As Long, ReadposA As Long, Len1 As Long, LastPos As Long, LenString As Long, MatrixLen As Long
    Dim NewAverage As Double, MultFactorx As Double, BB As Double, EE As Double
    Dim X As Integer, Y  As Integer, z  As Integer, A  As Integer, B As Integer
    Dim Boots() As String
    Dim MatrixByte() As Byte
    Dim LTree As Long, StartT As Long, NumWins As Long
    Dim RndNum2 As Long, RndNum As Long
    Dim NameLen As Integer, ttSeqNum() As Integer, Px() As Integer, xx1() As Integer, xx2() As Integer
    Dim WeightMod() As Long, Num1() As Long, Num2() As Long, Weight() As Long, Location() As Long, Ally() As Long, Alias() As Long
    Dim TMat2() As Double, Num() As Double, DistVal() As Double, Prod1() As Double, Prod2() As Double, Prod3() As Double, DEN() As Long
    Dim SHolder() As Byte
    Dim TreeOut As String
    Dim NodeOrder() As Long
    Dim NodeLen() As Double
    Dim DoneNode() As Long
    Dim TempNodeOrder() As Long
    Dim DstOut() As Integer, tSeqNum() As Integer
    Dim RootNode() As Long
    Dim TotMat2() As Double
    
    ReDim xx1(3)
    ReDim xx2(3)
    ReDim Prod1(Int(TOWinLen / 2))
    ReDim Prod2(Int(TOWinLen / 2))
    ReDim Prod3(Int(TOWinLen / 2))
    ReDim Alias(Int(TOWinLen / 2))
    ReDim Ally(Int(TOWinLen / 2))
    ReDim Location(Int(TOWinLen / 2))
    ReDim Px(ToNumSeqs, Int(TOWinLen / 2))
    ReDim TTempSeq2(Int(TOWinLen / 2), ToNumSeqs)
    ReDim DistVal(1)
    ReDim Num1(1)
    ReDim Num2(1)
    ReDim DEN(1)
    ReDim Num(1)
    ReDim WeightMod(0, Int(TOWinLen / 2))
    ReDim TotMat(ToNumSeqs, ToNumSeqs)
    ReDim SHolder((ToNumSeqs) * 40 * 2)
    ReDim Weight(0, Int(TOWinLen / 2))
    ReDim Boots(ToNumSeqs + 1)
    SeqLen = Len(StrainSeq(0))
    Dim NextWrite As Long
    If SeqLen > 32000 Then SeqLen = 32000
    Form1.ProgressBar1.Value = 5
    FF = 0
    CPos = 0
    Dim TempSeq2() As Integer
    ReDim TempSeq2(Int(TOWinLen / 2), ToNumSeqs)
    ReDim TotMat(ToNumSeqs, ToNumSeqs)
    Dim TS As String
    f2 = FreeFile
        On Error Resume Next
        Kill "intree"
        On Error GoTo 0
        Open "intree" For Binary Access Write As #f2
        
    NextWrite = 1
    eee = 0
    For X = 1 To TOPerms
        SSS = GetTickCount
        If Abs(SSS - eee) > 500 Then
            eee = SSS
            Form1.SSPanel1.Caption = "Simulating " & X & " of " & TOPerms & " Alignments"
            Form1.ProgressBar1.Value = 10 + (X / TOPerms) * 42
        End If
        BB = GetTickCount
        'Generate a simulated dataset using seqgen
        On Error Resume Next
            Kill "outfilex"
            Kill "out"
        On Error GoTo 0
        Call MakeSeqGenBat2(X)
        
        
        ShellAndClose "seqgen.bat", 0

        If AbortFlag = 1 Then
            Exit Sub
        End If

        'Load simulated datasets into a string array (ToSeq)
        FF2 = FreeFile
        Open "outfilex" For Input As #FF2

        If LOF(FF2) = 0 Then
            MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be SEQGEN.")
            AbortFlag = 1
            Close #FF2
            Close FF
            Exit Sub
        End If

        Line Input #FF2, Header

        If Len(Header) = LOF(FF2) Then 'some systems will not simply read lines - they go for the entire file
            Close #FF2
            Target = Chr$(10)
            LastPos = 1

            Do
                Pos = InStr(LastPos, Header, Target, vbBinaryCompare)

                If Pos > 0 Then
                    Header = left$(Header, Pos - 1) + Chr$(13) + right$(Header, Len(Header) - (Pos - 1))
                    LastPos = Pos + 3
                Else
                    Exit Do
                End If

            Loop

            Open "outfilex" For Output As #FF2
            Print #FF2, Header
            Close #FF2
            Open "outfilex" For Input As #FF2
            Line Input #FF2, Header
        End If

        For Y = 0 To ToNumSeqs
            Header = ""

            Do Until Len(Header) > 0
                Line Input #FF2, Header
            Loop
            z = Val(Trim$(Mid$(Header, 2, 5)))
            TOSeq(z) = right$(Header, Len(Header) - 10)
        Next 'Y

        Close #FF2
        
        'Name "outfilex" As "outfile" & X
        'convert toseq to integer array
        
        Dim TSeqSpaces() As Integer
        ReDim TSeqSpaces(Len(TOSeq(0)), ToNumSeqs)
        For z = 0 To ToNumSeqs
            Dummy = CopyString(Len(TOSeq(0)), TOTSeqNum(0, z), TOSeq(z), TSeqSpaces(0, z))
        Next 'Z
        
        'Do the windows

        If FF = 0 Then
            FF = FreeFile
            Open "infile" For Output As #FF
        End If

        StartPos = 1
        
        
        Do While StartPos + TOWinLen < SeqLen
            
            'Do the first half of the window.
            MakeSubAlign ToNumSeqs, Len(TOSeq(0)), StartPos + Int(TOWinLen / 4), 1, Int(TOWinLen / 2), TempSeq2(0, 0), TOTSeqNum(0, 0)
            DNADIST TOCoeffVar, TOTvTs, TOFreqFlag, TOModel, TOFreqA, TOFreqC, TOFreqG, TOFreqT, ToNumSeqs + 1, Int(TOWinLen / 2), TempSeq2(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), TotMat(0, 0)
            
            Call AverageMatrix(ToNumSeqs, NewAverage)

            If NewAverage = 0 Then
                MultFactorx = 1
            Else
                MultFactorx = MatAveragex / NewAverage
            End If

            Call MatrixMultiply(ToNumSeqs, MultFactorx)
            
            Call DistsToFile(FF, ToNumSeqs, TotMat())
            
            ReDim TotMat2(ToNumSeqs, ToNumSeqs)
            For z = 0 To ToNumSeqs - 1
                For Y = z + 1 To ToNumSeqs
                    'If TotMat(Z, Y) = 0 Then TotMat(Z, Y) = 0.00005
                    TotMat2(z, Y) = TotMat(z, Y)
                    TotMat2(Y, z) = TotMat(z, Y)
                Next Y
            Next z
            
            Pos = 0
            
            'This is where the writing to file thing used to go
            ReDim SHolder((ToNumSeqs + 1) * 40 * 2)
            'Do NJ tree
            ReDim ColTotals(ToNumSeqs)
            
            ReDim TreeArray(ToNumSeqs, ToNumSeqs)
            LTree = NEIGHBOUR(1, 2, TORndNum, 1, ToNumSeqs + 1, TotMat(0, 0), SHolder(0), ColTotals(0), TreeArray(0, 0))
            LTree = LTree + 2
            SHolder(LTree - 1) = 10
            SHolder(LTree) = 13
            Put #f2, NextWrite, SHolder
            NextWrite = NextWrite + LTree
            'Do second half of window
            MakeSubAlign ToNumSeqs, Len(TOSeq(0)), StartPos + Int(TOWinLen / 2) + Int(TOWinLen / 4), 1, Int(TOWinLen / 2), TempSeq2(0, 0), TOTSeqNum(0, 0)
            DNADIST TOCoeffVar, TOTvTs, TOFreqFlag, TOModel, TOFreqA, TOFreqC, TOFreqG, TOFreqT, ToNumSeqs + 1, Int(TOWinLen / 2), TempSeq2(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), TotMat(0, 0)
            
            Call AverageMatrix(ToNumSeqs, NewAverage)

            If NewAverage = 0 Then
                MultFactorx = 1
            Else
                MultFactorx = MatAveragex / NewAverage
            End If

            Call MatrixMultiply(ToNumSeqs, MultFactorx)
            
            Call DistsToFile(FF, ToNumSeqs, TotMat())
            
            
            ReDim TotMat2(ToNumSeqs, ToNumSeqs)
            For z = 0 To ToNumSeqs - 1
                For Y = z + 1 To ToNumSeqs
                    'If TotMat(Z, Y) = 0 Then TotMat(Z, Y) = 0.00005
                    TotMat2(z, Y) = TotMat(z, Y)
                    TotMat2(Y, z) = TotMat(z, Y)
                Next Y
            Next z
            
            Pos = 0
            
            'This is where the writing to file thing used to go
            ReDim SHolder((ToNumSeqs + 1) * 40 * 2)
            'Do NJ tree
            ReDim ColTotals(ToNumSeqs)
            
            ReDim TreeArray(ToNumSeqs, ToNumSeqs)
            LTree = NEIGHBOUR(1, 2, TORndNum, 1, ToNumSeqs + 1, TotMat(0, 0), SHolder(0), ColTotals(0), TreeArray(0, 0))
            LTree = LTree + 2
            SHolder(LTree - 1) = 10
            SHolder(LTree) = 13
            Put #f2, NextWrite, SHolder
            NextWrite = NextWrite + LTree
            CPos = CPos + 1
            StartPos = StartPos + TOStepSize
        Loop

        '3.6 seconds per perm
        
    Next 'X

    Close #FF
    Close #f2
    NumDatasets = (CPos - 1) * 2
    'Now we can sort out the optfiles
    
    'Fitch1 option file
    Open "optfilef1" For Output As #1
    Print #1, "g"
    Print #1, "j"
    Print #1, TORndNum
    Print #1, "1"
    Print #1, "p"
    Print #1, TOPower
    Print #1, "m"
    Print #1, NumDatasets
    Print #1, "3"
    Print #1, "3"
    Print #1, "2"
    Print #1, "y"
    Close #1
    'Fitch2 option file - used with NJ generated tree
    'and SS genration with "user" defined tree
    Open "optfilef2" For Output As #1
    Print #1, "u"
    Print #1, "p"
    Print #1, TOPower
    Print #1, "m"
    Print #1, NumDatasets
    Print #1, "3"
    Print #1, "2"
    Print #1, "3"
    Print #1, "y"
    Close #1
    Open "optfilen1" For Output As #1
    Print #1, "j"
    Print #1, TORndNum
    Print #1, "2"
    Print #1, "m"
    Print #1, NumDatasets
    Print #1, "2"
    Print #1, "3"
    Print #1, "y"
    Close #1
    
    If AbortFlag = 1 Then
        Exit Sub
    End If

    Form1.ProgressBar1.Value = 58
    '4.7 seconds
    StartPos = 1
    CPos = 0
    

    

    
    StartPos = 1
    Form1.SSPanel1.Caption = "Calculating SS Scores"
    CPos = 0
    Form1.ProgressBar1.Value = 63
    '3.7 seconds
    BB = GetTickCount
    If X = 12345 Then
        If TOTreeType = 1 Or TOPFlag = 1 Then  'do NJ if it is either set or you're working with permutations.
            
            
            ShellAndClose "neighbor1.bat", 0
            
            If AbortFlag = 1 Then
                Exit Sub
            End If
    
            
            FileCopy "infile", "infilebak"
            
            Kill "intree"
            Name "treefile" As "intree"
            
            'passing NJ calculated tree to fitch within the batch
            'infile requires some moving of trees and distance matrices
            'about
            
            On Error Resume Next
            Kill "outfile"
            Kill "outtree"
            On Error GoTo 0
            SS = GetTickCount
            ExpectFL = NumDatasets * 30: BatIndex = 9
            StartProgress = Form1.ProgressBar1.Value
            EndProgress = 90
            ShellAndClose "fitch2.bat", 0

            If AbortFlag = 1 Then
                Exit Sub
            End If
    
            Kill "infile"
            Name "infilebak" As "infile"
            '22 seconds
        Else
            On Error Resume Next
            Kill "outfile"
            Kill "outtree"
            On Error GoTo 0
            ExpectFL = NumDatasets * 30: BatIndex = 9
            StartProgress = Form1.ProgressBar1.Value
            EndProgress = 90
            ShellAndClose "fitch1.bat", 0
    
            If AbortFlag = 1 Then
                Exit Sub
            End If
    
        End If
    Else
        FileCopy "infile", "infilebak"
            
           ' Kill "intree"
           ' Name "treefile" As "intree"
            
            'passing NJ calculated tree to fitch within the batch
            'infile requires some moving of trees and distance matrices
            'about
            
            On Error Resume Next
            Kill "outfile"
            Kill "outtree"
            On Error GoTo 0
            SS = GetTickCount
            ExpectFL = NumDatasets * 30: BatIndex = 9
            StartProgress = Form1.ProgressBar1.Value
            EndProgress = 90
            ShellAndClose "fitch2.bat", 0

            If AbortFlag = 1 Then
                Exit Sub
            End If
    
            Kill "infile"
            Name "infilebak" As "infile"
    End If
    EE = GetTickCount
    TT = EE - BB
    Open "outfilex" For Binary Access Read As #1

    If LOF(1) = 0 Then
        MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component FITCH.")
        AbortFlag = 1
        Close #1
        Exit Sub
    End If

    GetStringA = String$(LOF(1), " ")
    Get #1, 1, GetStringA
    Close #1
    Form1.ProgressBar1.Value = 90
    'get ss scores for best fit trees
    ReDim SSScore(3, NumDatasets + 2)
    Pos = 1
    CPos = 0

    Do
        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            'SSScore(0, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error Resume Next
            SSScore(0, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(0, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
        End If

        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            
            On Error Resume Next
            SSScore(2, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(2, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
            CPos = CPos + 1
        End If

    Loop

    Form1.SSPanel1.Caption = "Recalculating SS Scores"
    Open "infile" For Binary Access Read As #1

    If LOF(1) = 0 Then
        MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component FITCH.")
        AbortFlag = 1
        Close #1
        Exit Sub
    End If

    GetStringA = String$(LOF(1), " ")
    Get #1, 1, GetStringA
    Close #1
    
    'Work out distance matrix spacing
    Pos = InStr(1, GetStringA, "S0", vbBinaryCompare)
    HeaderLen = Pos - 1
    LastPos = Pos
    Pos = InStr(LastPos + 1, GetStringA, "S00", vbBinaryCompare)
    Len1 = Pos - LastPos
    LenString = Len(GetStringA)
    Pos = 1
    LastPos = -2
    ReadposA = 1
    'Write the infile for SS calculation with fitch
    Open "infile" For Output As #1
    
    NumDatasets = 0
    Do While (ReadposA + Len1 * 2) - 1 <= LenString
        'get distances
        GetStringC = Mid$(GetStringA, ReadposA, Len1)
        ReadposA = ReadposA + Len1
        GetStringE = Mid$(GetStringA, ReadposA, Len1)
        ReadposA = ReadposA + Len1
        
        'note:previous cycle's trees used
       
        Print #1, GetStringE
        
        Print #1, GetStringC
        
        NumDatasets = NumDatasets + 2
    Loop
    
    Close #1
    On Error Resume Next
        Kill "outfile"
        Kill "outtree"
        On Error GoTo 0
        SS = GetTickCount
    ExpectFL = NumDatasets * 30: BatIndex = 9
    StartProgress = Form1.ProgressBar1.Value
    EndProgress = 98
    ShellAndClose "fitch2.bat", 0
    EE = GetTickCount
    TT = EE - SS
    If AbortFlag = 1 Then
        Exit Sub
    End If

    Form1.ProgressBar1.Value = 98
    '4.6 seconds
    BB = GetTickCount
    Open "outfilex" For Binary Access Read As #1

    If LOF(1) = 0 Then
        MsgBox ("An error occured during the execution of TOPAL.  The problem appears to be the Phyip component FITCH.")
        AbortFlag = 1
        Close #1
        Exit Sub
    End If

    GetStringA = String$(LOF(1), " ")
    Get #1, 1, GetStringA
    Close #1
    'get ssscores
    Pos = 1
    CPos = 0

    Do
        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
           
            On Error Resume Next
            SSScore(1, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(1, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
        End If

        Pos = InStr(Pos, GetStringA, "Sum of squares =   ", vbBinaryCompare)

        If Pos = 0 Then
            Exit Do
        Else
            
            On Error Resume Next
            SSScore(3, CPos) = Val(Mid$(GetStringA, Pos + 21, 7))
            SSScore(3, CPos) = CDbl(Mid$(GetStringA, Pos + 21, 7))
            On Error GoTo 0
            Pos = Pos + 1
            CPos = CPos + 1
        End If

    Loop
    '0.2 seconds
    'CPos = NumDatasets
End Sub



Public Sub MakeSeqGenBat2(CurPerm)
'Write seq-gen batch file
Dim AF As Double, CF As Double, GF As Double, TF As Double
    FFX = FreeFile
    Open "seqgen.bat" For Output As #FFX
    'different cline parametes must be given for different models
    'onlt JC model options are given
    'Note that for JN and ML extra things must be either calculated
    'Shape parameter alpha distribution for JN and bpfrequencies for
    'ML).
    ToSeqGenCLine = "seq-gen"

    If TOModel = 0 Then
        ToSeqGenCLine = ToSeqGenCLine & " -mHKY -t0.5"
    ElseIf TOModel = 1 Then
        ToSeqGenCLine = ToSeqGenCLine & " -mHKY -t" & TOTvTs
    ElseIf TOModel = 2 Then
        ToSeqGenCLine = ToSeqGenCLine & " -mHKY -t" & TOTvTs
    Else
        ToSeqGenCLine = ToSeqGenCLine & " -mF84 -t" & TOTvTs

        If TOFreqFlag = 0 Then

            Call CalcBPFreqs(AF, CF, GF, TF)

            ToSeqGenCLine = ToSeqGenCLine & " -f " & AF & " " & CF & " " & GF & " " & TF
        Else
            ToSeqGenCLine = ToSeqGenCLine & " -f " & TOFreqA & " " & TOFreqC & " " & TOFreqG & " " & TOFreqT
        End If

    End If
    ToSeqGenCLine = ToSeqGenCLine & " -z" & CStr(CurPerm)
    
    If Len(TOSeq(0)) < 32000 Then
        ToSeqGenCLine = ToSeqGenCLine & " -l" & Len(TOSeq(0)) & " -n1 <simtree > out"
    Else
        ToSeqGenCLine = ToSeqGenCLine & " -l32000 -n1 <simtree > out"
    End If
    
    
    
    Print #FFX, ToSeqGenCLine
    Print #FFX, "del outfilex"
    Print #FFX, "rename out outfilex"
    Close #FFX
    
End Sub
Public Sub LXoverA()
 'This subroutine executes LARD for automated screens - it only uses 2 bp scans with a heuristic search
    Dim PVal As Double, LRDRegion As Byte, BPos As Long, EPos As Long, LRDWin
    LRDRegion = 2
    LRDWin = 0
    Dim TSeq(2) As String, X As Long
    For X = 1 To Len(StrainSeq(0))
        
        TSeq(0) = TSeq(0) + Chr(SeqNum(X, Seq1) - 1)
        TSeq(1) = TSeq(1) + Chr(SeqNum(X, Seq2) - 1)
        TSeq(2) = TSeq(2) + Chr(SeqNum(X, Seq3) - 1)
        
    Next X
    
    Dim PID As Long, ProcessID As Long, LastLen As Long, SCount As Long, FLen As Long, NewSurface As Long, Pict As Long, CurPos As Long, Count As Long, RetVal As Long, NewPos As Long
    Dim PntAPI As POINTAPI
    Dim OldDir As String, LARDCLine As String, TitleTmp As String

    Const STILL_ACTIVE = &H103
    ReDim LXPos(Len(StrainSeq(0)))
    

    OldDir = CurDir
    ChDir App.Path
    ChDrive App.Path
    
    
    'Writes the alignment file for LARD
    Open "lardin" For Output As #1
    
    
    
    'LRDWin = 1
    'LRDWinLen = 200
    
    
    If LRDWin = 0 Then
        
        NumWins = 1
        
        Print #1, " 3  " + CStr(Len(StrainSeq(0))) ' 200"
        Print #1, "s1"
        Print #1, TSeq(0)
        Print #1, "s2"
        Print #1, TSeq(1)
        Print #1, "s3"
        Print #1, TSeq(2)
        
    Else
        
        'LRDWinLen = 200
        NumWins = (Len(StrainSeq(0)) - LRDWinLen) / LRDStep
        
        For X = 0 To NumWins - 1
            Print #1, " 3  " + CStr(LRDWinLen)
            z = X * LRDStep + 1
            Print #1, "s1"
            Print #1, Mid$(TSeq(0), z, LRDWinLen)
            Print #1, "s2"
            Print #1, Mid$(TSeq(1), z, LRDWinLen)
            Print #1, "s3"
            Print #1, Mid$(TSeq(2), z, LRDWinLen)
            Print #1, ""
        Next X
    End If
    Close #1
    'Exit Sub
    'Creates a fake outfile and then kills it
    Open "likelihood.surface" For Binary As 1
    Put #1, 1, ""
    Close #1
    On Error Resume Next
    Kill "likelihood.surface"
    On Error GoTo 0
    'Creates a fake stackdump file and then kills it
    Open "LARD.EXE.stackdump" For Output As #1
    Close #1
    Kill "LARD.EXE.stackdump"
    'Constructs the LARD command line
    LARDCLine = "lard.exe"
    
    If LRDWin = 0 Then
        LARDCLine = LARDCLine + " -r" + Trim(Str(CInt(LRDRegion)))
    Else
        LARDCLine = LARDCLine + " -r1"
    End If
    If LRDModel = 0 Or LRDModel = 1 Then

        If LRDModel = 0 Then
            LARDCLine = LARDCLine & " -mHKY"
        Else
            LARDCLine = LARDCLine & " -mF84"
        End If

        LARDCLine = LARDCLine & " -t" & LRDTvRat

        If LRDModel = 0 Or (LRDModel = 1 And LRDBaseFreqFlag = 1) Then

            If LRDAFreq <> LRDCFreq Or LRDAFreq <> LRDGFreq Or LRDAFreq <> LRDTFreq Then
                LARDCLine = LARDCLine & " -f" & LRDAFreq & " " & LRDCFreq & " " & LRDGFreq & " " & LRDTFreq
            End If

        End If

    Else
        LARDCLine = LARDCLine & " -mREV"

        If LRDAFreq <> LRDCFreq Or LRDAFreq <> LRDGFreq Or LRDAFreq <> LRDTFreq Then
            LARDCLine = LARDCLine & " -f" & LRDAFreq & " " & LRDCFreq & " " & LRDGFreq & " " & LRDTFreq
        End If

        If LRDACCoeff <> 1 Or LRDAGCoeff <> 1 Or LRDATCoeff <> 1 Or LRDCGCoeff <> 1 Or LRDCTCoeff <> 1 Then
            LARDCLine = LARDCLine & " -t" & LRDACCoeff & " " & LRDAGCoeff & " " & LRDATCoeff & " " & LRDCGCoeff & " " & LRDCTCoeff
        End If

    End If
    lrdgd = 0
    
    If LRDWin = 1 Then
        LARDCLine = LARDCLine + " -n" + CStr(Int(NumWins))
        LARDCLine = LARDCLine + " -u" + CStr(CLng(LRDWinLen / 2))
    Else
        If lrdgd = 0 Then
            LARDCLine = LARDCLine & " -e -s" & LRDStep
        Else
            LARDCLine = LARDCLine & " -d -s" & LRDStep
        End If
        
        LARDCLine = LARDCLine & " -z" + Trim(Str(CInt(LRDStep * 2)))
    End If
    
    
   

    If LRDCodon1 <> LRDCodon2 Or LRDCodon1 <> LRDCodon3 Then
        LARDCLine = LARDCLine & " -c" & LRDCodon1 & " " & LRDCodon2 & " " & LRDCodon3
    End If

    If LRDCategs > 0 Then
        LARDCLine = LARDCLine & " -g" & LRDCategs
        LARDCLine = LARDCLine & " -a" & LRDShape
    End If
    LARDCLine = LARDCLine & " -vp"
    LARDCLine = LARDCLine & " <lardin"
    'Creates a batch file so that lard can be executed with a command line
    
    
    'LARDCLine = "cmd.exe /c " + LARDCLine + " >c:\test.txt"
    'LARDCLine = "cmd.exe /c ping /? >c:\test.txt"
    X = X
    Open "lard.bat" For Output As #1
    
    Do
        Pos = InStr(1, LARDCLine, ",")
        If Pos = 0 Then Exit Do
        Mid$(LARDCLine, Pos, 1) = "."
    Loop
    Print #1, LARDCLine
    Close #1
    'Checks to see if a previously aborted run of lard was made and terminates it if necessary

    
    'ProcessID = Shell("lard.bat ", 0)
    'Do
    
    Dim OP As String
    'ProcessID = Shell("lard.bat", 1)
    'OP = GetCommandOutput(LARDCLine, True, True)
    'SS = GetTickCount
    
    Form1.Command25.Enabled = True
    
    Dim EC As Long
    EC = CInt(Len(StrainSeq(0)) / LRDStep) + 2
    'pvalcalc = 2*lik ratio chisquare with 4 df per breakpoint
    If LRDRegion = 2 Then
        
        ExpectFL = EC * EC * 22
            
    End If
    
    If LRDWin = 1 Then
        OP = GetCommandOutput("lard.bat", 0, True, False)
    Else
        OP = GetCommandOutput("lard.bat", 0, True, True)
    End If
    Open "screen.out" For Output As #1
    Print #1, OP
    Close #1
    'check to see if this beats the p-val cutoff
    If MCFlag = 0 Then
        MC = MCCorrection * (Len(StrainSeq(0)) / LRDStep)
    Else
        MC = (Len(StrainSeq(0)) / LRDStep)
    End If
    
    
    
    Pos = InStr(1, OP, "LR=", vbBinaryCompare)
    Pos = Pos + 3
    Pos2 = InStr(Pos + 1, OP, Chr(13), vbBinaryCompare)
    If Pos > 22 Then
        MaxL = Mid$(OP, Pos, Pos2 - Pos)
    
    End If
    If MaxL > 745 Then MaxL = 745
    PVal = chi2(MaxL * 2, 4) * MC
    
    
    If PVal > LowestProb Then Exit Sub
    
    If LRDWin = 0 Then
    
        'there is no likelihood surface file with heuristic screen so don't bother with that
        'read likelihoods from the OP output
        Pos = InStr(1, OP, "Recombination region (to the left of):", vbBinaryCompare)
        If Pos = 0 Then Exit Sub
        LastPos = Pos + 1
        
        Pos = InStr(LastPos, OP, "->", vbBinaryCompare)
        XX = Mid$(OP, LastPos + 38, Pos - LastPos - 38)
        BPos = CLng(Mid$(OP, LastPos + 38, Pos - LastPos - 38))
        LastPos = Pos + 2
        Pos = InStr(LastPos, OP, Chr(9), vbBinaryCompare)
        EPos = CLng(Mid$(OP, LastPos, Pos - LastPos))
        X = X
    Else
        'work out bpos and epos however I can - possibly split the alignment at maxpeak and do single scan for next maxpeak.
        'read positions from op
        ReDim LXPos(0, NumWins)
        ReDim LSurface(NumWins)
       
        X = X
        Pos = InStr(1, OP, "cross", vbBinaryCompare)
        If Pos > 0 Then
            LastPos = Pos
            
            For X = 0 To Int(NumWins) - 1
                For Y = 0 To 2
                    Pos = InStr(LastPos + 1, OP, Chr(9), vbBinaryCompare)
                    LastPos = Pos
                    XX = Mid$(OP, Pos, 20)
                    X = X
                Next Y
                TPos = Pos
                Pos = InStr(LastPos + 1, OP, Chr(9), vbBinaryCompare)
                LastPos = Pos
                
                LSurface(X + 1) = CDbl(Mid$(OP, TPos + 1, Pos - TPos - 1))
                LXPos(0, X + 1) = (LRDStep * X) + LRDWinLen / 2
            Next X
        Else
            'handle this sometime
        End If
    End If
        
    
    Dim ScoresX() As Byte
    ReDim ScoresX(Len(StrainSeq(0)), 2)
    LenXoverSeq = BSSubSeq(Len(StrainSeq(0)), Seq1, Seq2, Seq3, SeqNum(0, 0), XPosDiff(0), XDiffpos(0), ScoresX(0, 0))
    
    Call ProcessEvent(7, PVal, BPos, EPos, XPosDiff(), XDiffpos(), 0, 0, 10, LenXoverSeq)
    
    Call ProcessEvent(7, PVal, EPos, BPos, XPosDiff(), XDiffpos(), 0, 0, 10, LenXoverSeq)
    'check to see whether the best even is significant and if it is add it to xoverlist with a prgflag=6
    
        
    

End Sub
Public Sub AddToBak(SuperEventList() As Long, ReplaceE, WinE, oEventNumber, S2TraceBack() As Long, oSuperEventList() As Long, BakXOList() As XOverDefine, BakCurXOver() As Integer, XOverList() As XOverDefine, CurrentXOver() As Integer)
If ReplaceE > 1 Then
    X = X '1471,1472,1471,1469:1498,1471,1472,1471,1469
End If
For X = 0 To Nextno
    For Y = 1 To CurrentXOver(X)
        If WinE = SuperEventList(XOverList(X, Y).Eventnumber) Then
                    
            BakCurXOver(S2TraceBack(X)) = BakCurXOver(S2TraceBack(X)) + 1
            If BakCurXOver(S2TraceBack(X)) > UBound(BakXOList, 2) Then
                ReDim Preserve BakXOList(PermNextNo, BakCurXOver(S2TraceBack(X)) + 100)
            End If
            
            BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))) = XOverList(X, Y)
            BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).Daughter = S2TraceBack(XOverList(X, Y).Daughter)
            BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).MajorP = S2TraceBack(XOverList(X, Y).MajorP)
            BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).MinorP = S2TraceBack(XOverList(X, Y).MinorP)
            XX = XOverList(X, Y).Beginning
            BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).Eventnumber = oEventNumber + XOverList(X, Y).Eventnumber
            
            If BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).Eventnumber > UBound(oSuperEventList, 1) Then
                ReDim Preserve oSuperEventList(BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).Eventnumber + 100)
            End If
            oSuperEventList(BakXOList(S2TraceBack(X), BakCurXOver(S2TraceBack(X))).Eventnumber) = ReplaceE
            
        End If
    Next Y
Next X
'XX = BakXOlist(1, 6).Eventnumber

End Sub

Public Sub GetAge(EN, Daught() As Long, G, d, P1, P2, FMat() As Double, SMat() As Double, PermDiffs() As Double, PermValid() As Double, AgeEvent() As Double)
Dim X As Long, Y As Long
Dim Age As Double, Addj(1) As Double, ID() As Double, Tot(2) As Double, TParDist As Double, ParDist As Double
Dim MaxDist() As Double
If SEventNumber = 1 Then
    ReDim AgeEvent(1, 10)
Else
    If UBound(AgeEvent, 1) < SEventNumber Then
        ReDim Preserve AgeEvent(1, SEventNumber + 10)
    End If
End If
ReDim ID(1, 2)
ID(0, 0) = FMat(d, P1): ID(0, 1) = FMat(d, P2): ID(0, 2) = FMat(P1, P2)
ID(1, 0) = SMat(d, P1): ID(1, 1) = SMat(d, P2): ID(1, 2) = SMat(P1, P2)
'Work out the distance modifyer for each region
Tot(0) = 0: Tot(1) = 0: Tot(2) = 0
For X = 0 To Nextno
    For Y = X + 1 To Nextno
        'If DontUse(X) = 0 And DontUse(Y) = 0 Then
            If FMat(X, Y) < 3 Then
                Tot(0) = Tot(0) + FMat(X, Y)
                Tot(1) = Tot(1) + SMat(X, Y)
                TParDist = 1 - PermDiffs(X, Y) / PermValid(X, Y)
                If TParDist > 0.25 Then
                    TParDist = (4# * TParDist - 1#) / 3#
                    TParDist = Log(TParDist)
                    TParDist = -0.75 * TParDist
                   
                Else
                    TParDist = 1
                End If
                Tot(2) = Tot(2) + TParDist
            End If
        'End If
    Next Y
Next X

If Tot(0) > 0 And Tot(2) > 0 Then '42.6091706,18.498428
    Addj(0) = Tot(0) / Tot(2)
Else
    Addj(0) = 1
End If
If Tot(1) > 0 And Tot(2) > 0 Then '42.6091706,18.498428
    Addj(1) = Tot(1) / Tot(2)
Else
    Addj(1) = 1
End If
'Modify the recombinant region dists
For X = 0 To 2
    ID(0, X) = ID(0, X) / Addj(0)
    ID(1, X) = ID(1, X) / Addj(1)
   
Next X

ParDist = 0: Age = 1000
For X = 0 To 1
    For Y = 0 To 2
        If Age > ID(X, Y) Then Age = ID(X, Y)
        If ParDist < ID(X, Y) Then ParDist = ID(X, Y)
        
    Next Y
Next X



ParDist = ParDist - Age



If ParDist < 0 Then ParDist = 0
'AgeEvent(1, G) = ParDist
AgeEvent(0, EN) = Age
'now work out the minimum age of each event
ReDim MaxDist(1)
For Y = 0 To Nextno
    If Daught(G, Y) > 0 Then
         If FMat(d, Y) > MaxDist(0) Then MaxDist(0) = FMat(d, Y)
         If SMat(d, Y) > MaxDist(1) Then MaxDist(1) = SMat(d, Y)
    End If
Next Y
MaxDist(0) = MaxDist(0) / Addj(0)
MaxDist(1) = MaxDist(1) / Addj(1)
If MaxDist(0) < MaxDist(1) Then
    AgeEvent(1, EN) = MaxDist(0)
Else
    AgeEvent(1, EN) = MaxDist(1)
End If

If AgeEvent(0, EN) < AgeEvent(1, EN) Then
    AgeEvent(1, EN) = AgeEvent(0, EN)
Else
    If AgeEvent(1, EN) = 0 Then
        AgeEvent(1, EN) = AgeEvent(0, EN) / 2
    Else
         AgeEvent(1, EN) = AgeEvent(1, EN) + (AgeEvent(0, EN) - AgeEvent(1, EN)) / 2
         
    End If
End If
X = X
End Sub

Public Sub EraseEvidence(Nextno, CE, BakXOList() As XOverDefine, BakCurXOver() As Integer, SuperEventList() As Long)
Dim X As Long, Y As Long
For X = 0 To Nextno
    Y = 1
    Do While Y <= BakCurXOver(X)
        If SuperEventList(BakXOList(X, Y).Eventnumber) = CE Then
            If Y < BakCurXOver(X) Then
                BakXOList(X, Y) = BakXOList(X, BakCurXOver(X))
            End If
            BakCurXOver(X) = BakCurXOver(X) - 1
        
        End If
        Y = Y + 1
    Loop
Next X
End Sub
Public Sub MakeBestEvent()

Dim BestP() As Double, AProg() As Byte
ReDim AProg(AddNum * 2)
'Eventnumber = 1000
If DoScans(0, 0) = 1 Then AProg(0) = 1: AProg(0 + AddNum) = 1
If DoScans(0, 1) = 1 Then AProg(1) = 1: AProg(1 + AddNum) = 1
If DoScans(0, 2) = 1 Then AProg(2) = 1: AProg(2 + AddNum) = 1
If DoScans(0, 3) = 1 Then AProg(3) = 1: AProg(3 + AddNum) = 1
If DoScans(0, 4) = 1 Then AProg(4) = 1: AProg(4 + AddNum) = 1
If DoScans(0, 5) = 1 Then AProg(5) = 1: AProg(5 + AddNum) = 1

ReDim BestEvent(SEventNumber, 1), BestP(SEventNumber)
ReDim Confirm(Eventnumber + 1, AddNum - 1), ConfirmP(Eventnumber + 1, AddNum - 1), ConfirmMi(Eventnumber + 1, AddNum - 1), ConfirmPMi(Eventnumber + 1, AddNum - 1), ConfirmMa(Eventnumber + 1, AddNum - 1), ConfirmPMa(Eventnumber + 1, AddNum - 1)

For X = 0 To Nextno
    For Y = 1 To CurrentXOver(X)
    
        pf = XOverList(X, Y).ProgramFlag
        If pf <= AddNum - 1 Then
            If XOverList(X, Y).Probability > 0 And (XOverList(X, Y).Probability < BestP(SuperEventList(XOverList(X, Y).Eventnumber)) Or BestP(SuperEventList(XOverList(X, Y).Eventnumber)) = 0) Then
                If AProg(XOverList(X, Y).ProgramFlag) = 1 Then
                    BestP(SuperEventList(XOverList(X, Y).Eventnumber)) = XOverList(X, Y).Probability
                    
                End If
                
            End If
            If XOverList(X, Y).DHolder < 0 Then
               
                BestEvent(SuperEventList(XOverList(X, Y).Eventnumber), 0) = X
                BestEvent(SuperEventList(XOverList(X, Y).Eventnumber), 1) = Y
            End If
            If pf > AddNum - 1 Then pf = pf - AddNum
            'XX = UBound(Confirm, 1)
            Confirm(SuperEventList(XOverList(X, Y).Eventnumber), pf) = Confirm(SuperEventList(XOverList(X, Y).Eventnumber), pf) + 1
            ConfirmP(SuperEventList(XOverList(X, Y).Eventnumber), pf) = ConfirmP(SuperEventList(XOverList(X, Y).Eventnumber), pf) + -Log10(XOverList(X, Y).Probability)
        End If
    Next Y
Next X

    
    
 

 For X = 1 To SEventNumber
     BP = 10000
     If BestEvent(X, 0) = 0 And BestEvent(X, 1) = 0 Then
         'make the nextbest the best
         For z = 0 To Nextno
             For Y = 1 To CurrentXOver(z)
                 If SuperEventList(XOverList(z, Y).Eventnumber) = X Then
                     If AProg(XOverList(z, Y).ProgramFlag) = 1 Then
                         If XOverList(z, Y).Probability > 0 And XOverList(z, Y).Probability < BP Then
                             BestEvent(X, 0) = z: BestEvent(X, 1) = Y
                             BP = XOverList(z, Y).Probability
                             
                         End If
                     End If
                 End If
             Next Y
         Next z
         
         If BestEvent(X, 0) <> 0 Or BestEvent(X, 1) <> 0 Then
             
             If XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder > 0 Then
                 XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder = -XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder
             ElseIf XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder = 0 Then
                 XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder = -0.00001
             End If
         End If
     
     End If
 
 Next X
 For X = 1 To SEventNumber
     BP = 10000
     If BestEvent(X, 0) = 0 And BestEvent(X, 1) = 0 Then
         'make the nextbest the best
         For z = 0 To Nextno
             For Y = 1 To CurrentXOver(z)
                 If SuperEventList(XOverList(z, Y).Eventnumber) = X Then
                    
                     If XOverList(z, Y).Probability > 0 And XOverList(z, Y).Probability < BP Then
                         BestEvent(X, 0) = z: BestEvent(X, 1) = Y
                         BP = XOverList(z, Y).Probability
                        
                     End If
                 End If
             Next Y
         Next z
         
     End If
     If BestEvent(X, 0) <> 0 Or BestEvent(X, 1) <> 0 Then
         If XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder > 0 Then
             XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder = -XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder
         ElseIf XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder = 0 Then
             XOverList(BestEvent(X, 0), BestEvent(X, 1)).DHolder = -0.00001
         
         End If
     End If
 Next X
 For X = 0 To Nextno
     For Y = 1 To CurrentXOver(X)
         pf = XOverList(X, Y).ProgramFlag
         If pf > AddNum - 1 Then
             If XOverList(X, Y).Probability > 0 And (XOverList(X, Y).Probability < BestP(SuperEventList(XOverList(X, Y).Eventnumber)) Or BestP(SuperEventList(XOverList(X, Y).Eventnumber)) = 0) Then
                 If AProg(XOverList(X, Y).ProgramFlag) = 1 Then
                     If BestP(SuperEventList(XOverList(X, Y).Eventnumber)) = 0 Then
                         BestP(SuperEventList(XOverList(X, Y).Eventnumber)) = XOverList(X, Y).Probability
                        ' BestEvent(SuperEventlist(XOverList(X, Y).Eventnumber), 0) = X
                         'BestEvent(SuperEventlist(XOverList(X, Y).Eventnumber), 1) = Y
                     Else
                         If XOverList(BestEvent(SuperEventList(XOverList(X, Y).Eventnumber), 0), BestEvent(SuperEventList(XOverList(X, Y).Eventnumber), 1)).ProgramFlag > AddNum - 1 Then
                             BestP(SuperEventList(XOverList(X, Y).Eventnumber)) = XOverList(X, Y).Probability
                           '  BestEvent(SuperEventlist(XOverList(X, Y).Eventnumber), 0) = X
                           '  BestEvent(SuperEventlist(XOverList(X, Y).Eventnumber), 1) = Y
                         End If
                     End If
                 End If
             End If
             
             If pf > AddNum - 1 Then pf = pf - AddNum
             Confirm(SuperEventList(XOverList(X, Y).Eventnumber), pf) = Confirm(SuperEventList(XOverList(X, Y).Eventnumber), pf) + 1
             ConfirmP(SuperEventList(XOverList(X, Y).Eventnumber), pf) = ConfirmP(SuperEventList(XOverList(X, Y).Eventnumber), pf) + -Log10(XOverList(X, Y).Probability)
         End If
     Next Y
 Next X
 
 For X = 0 To PermNextNo
     For Y = 1 To BCurrentXoverMi(X)
         pf = BestXOListMi(X, Y).ProgramFlag
         If pf > AddNum - 1 Then pf = pf - AddNum
         ConfirmMi(SuperEventList(BestXOListMi(X, Y).Eventnumber), pf) = ConfirmMi(SuperEventList(BestXOListMi(X, Y).Eventnumber), pf) + 1
         ConfirmPMi(SuperEventList(BestXOListMi(X, Y).Eventnumber), pf) = ConfirmPMi(SuperEventList(BestXOListMi(X, Y).Eventnumber), pf) + -Log10(BestXOListMi(X, Y).Probability)
     Next Y
 Next X
 For X = 0 To PermNextNo
     For Y = 1 To BCurrentXoverMa(X)
         pf = BestXOListMa(X, Y).ProgramFlag
         If pf > AddNum - 1 Then pf = pf - AddNum
         ConfirmMa(SuperEventList(BestXOListMa(X, Y).Eventnumber), pf) = ConfirmMa(SuperEventList(BestXOListMa(X, Y).Eventnumber), pf) + 1
         ConfirmPMa(SuperEventList(BestXOListMa(X, Y).Eventnumber), pf) = ConfirmPMa(SuperEventList(BestXOListMa(X, Y).Eventnumber), pf) + -Log10(BestXOListMa(X, Y).Probability)
     Next Y
 Next X
End Sub
Public Sub Scan2()
'Exit Sub
Dim oSelGrpFlag As Byte, oGrpMask() As Byte, EN As Long, AgeOrder() As Long, StartNextno As Long, Trace() As Long, NewSeqnum() As Integer


XX = Eventnumber
StartNextno = Nextno
oSelGrpFlag = SelGrpFlag
SelGrpFlag = 1
ReDim S2TraceBack(Nextno), MissingData(Len(StrainSeq(0)), Nextno), oGrpMask(Nextno)
For X = 0 To Nextno
    oGrpMask(X) = GrpMaskSeq(X)
    S2TraceBack(X) = X
Next X

ReDim Preserve TreeTrace(Nextno + 100), CurrentXOver(Nextno + 100), Daught(SEventNumber, Nextno + 100), StraiName(Nextno + 100), SeqNum(Len(StrainSeq(0)), Nextno + 100), MissingData(Len(StrainSeq(0)), Nextno + 100)

ReDim NewSeqnum(Len(StrainSeq(0)), Nextno)



If X = 12345 Then
    For X = 0 To StepNo
        
        If Steps(0, X) = 1 Then 'ie make a sequence
            
            
            Nextno = Nextno + 1
            UB = UBound(MissingData, 2)
            If Nextno > UB Then
                ReDim Preserve TreeTrace(Nextno + 100), CurrentXOver(Nextno + 100), Daught(SEventNumber, Nextno + 100), StraiName(Nextno + 100), SeqNum(Len(StrainSeq(0)), Nextno + 100), MissingData(Len(StrainSeq(0)), Nextno + 100)
            End If
                            
            StraiName(Nextno) = StraiName(TreeTrace(Steps(1, X))) 'x97704(6)
            Daught(EN, Nextno) = Daught(EN, Steps(1, X))
                            
            TreeTrace(Nextno) = TreeTrace(TreeTrace(Steps(1, X)))
                            
            Dummy = ModSeqNumD(Nextno, Len(StrainSeq(0)), Steps(1, X), Steps(2, X), Steps(3, X), SeqNum(0, 0), MissingData(0, 0))
                            
                            
            If Steps(2, X) < Steps(3, X) Then
                For Y = Steps(2, X) To Steps(3, X)
                    NewSeqnum(Y, TreeTrace(Steps(1, X))) = NewSeqnum(Y, TreeTrace(Steps(1, X))) + 1
                Next Y
            Else
                For Y = Steps(2, X) To Len(StrainSeq(0))
                    NewSeqnum(Y, TreeTrace(Steps(1, X))) = NewSeqnum(Y, TreeTrace(Steps(1, X))) + 1
                Next Y
                For Y = 1 To Steps(3, X)
                    NewSeqnum(Y, TreeTrace(Steps(1, X))) = NewSeqnum(Y, TreeTrace(Steps(1, X))) + 1
                Next Y
            End If
                            
        ElseIf Steps(0, X) = 2 Then  'delete some positions in sequence
                            
            Dummy = ModSeqNumE(Nextno, Len(StrainSeq(0)), Steps(1, X), Steps(2, X), Steps(3, X), SeqNum(0, 0), MissingData(0, 0))
                            
        ElseIf Steps(0, X) = 3 Then 'delete a sequence
            If Steps(1, X) < Nextno Then
                                
                Dummy = ReplaceSeq(Nextno, Len(StrainSeq(0)), Steps(1, X), Nextno, SeqNum(0, 0), MissingData(0, 0))
                Daught(EN, Steps(1, X)) = Daught(EN, Nextno)
                TreeTrace(Steps(1, X)) = TreeTrace(TreeTrace(Nextno))
                StraiName(Steps(1, X)) = StraiName(Nextno)
            End If
            Daught(EN, Nextno) = 0
            Nextno = Nextno - 1
        ElseIf Steps(0, X) = 4 Then 'replace a sequence
                            
            Dummy = ReplaceSeq(Nextno, Len(StrainSeq(0)), Steps(1, X), Nextno, SeqNum(0, 0), MissingData(0, 0))
            Daught(EN, Steps(1, X)) = Daught(EN, Nextno)
            StraiName(Steps(1, X)) = StraiName(Nextno)
                            
            TreeTrace(Steps(1, X)) = TreeTrace(TreeTrace(Nextno))
                            
            Nextno = Nextno - 1
        End If
    
    Next X
End If
Call UnModSeqNum(0)
Call UnModNextno

Dim TempXOList() As XOverDefine, DoneSeq() As Byte, CE As Long, AdjDst As Double, TParDist As Double, PairDist() As Double, PairValid As Double, PairDiff As Double, DateMat() As Double, SConvert As Double, ST As Long
Dim UseAll As Byte
ReDim TempXOList(Nextno, 0)
ReDim DoneSeq(Nextno, 1)



Dim MaxDim As Long, BakXOList() As XOverDefine, BakCurXOver() As Integer, BakXOlistMi() As XOverDefine, BakCurXOverMi() As Integer, BakXOlistMa() As XOverDefine, BakCurXOverMa() As Integer
ReDim Preserve Daught(SEventNumber, Nextno + 100)
CE = 1
Dim FSEN As Long
FSEN = SEventNumber


    Call UnModNextno
    Call UnModSeqNum(0)
     ReDim Preserve S2TraceBack(Nextno + 100)
    
    ReDim AgeEvent(1, SEventNumber), AgeOrder(SEventNumber)
    
    
    
    If X = X Then
        UseAll = 0
        Dim Size As Long, d As Long, P1 As Long, P2 As Long, Win As Long
        
        ENumb = SEventNumber
        
        Dim Excl() As Byte, Enu As Long, NC As Long
        ReDim Excl(SEventNumber), pCount(SEventNumber, AddNum)
            
        Dim BPV() As Double
        ReDim BPV(SEventNumber, AddNum)
            
        For X = 1 To SEventNumber
            Excl(X) = 1
        Next X
            
            
        zz = 0
                
    
        Dim Age As Double, RS1 As Long, LS1 As Long, RS2 As Long, LS2 As Long, NS(1) As Long, ParDist As Double
        xSeq1 = Seq1
        xSeq2 = Seq2
        xSeq3 = Seq3
        xrelx = RelX
        xrely = RelY
        xnjflag = NJFlag
        NJFlag = 0
        
        
        Call UnModSeqNum(0)
        
        For G = 1 To ENumb
            
            If BestEvent(G, 0) > 0 Or BestEvent(G, 1) > 0 Then
                
                CNum = 0
                sen = X 'SuperEventlist(XOverList(BestEvent(X, 0), BestEvent(X, 1)).Eventnumber)
                 
                    
                d = XOverList(BestEvent(G, 0), BestEvent(G, 1)).Daughter
                P1 = XOverList(BestEvent(G, 0), BestEvent(G, 1)).MajorP
                P2 = XOverList(BestEvent(G, 0), BestEvent(G, 1)).MinorP
                ST = XOverList(BestEvent(G, 0), BestEvent(G, 1)).Beginning
                EN = XOverList(BestEvent(G, 0), BestEvent(G, 1)).Ending
                
                If X = X Then
                    If X = X Then
                    
                        Seq1 = P1
                        Seq2 = P2
                        Seq3 = d
                        RelX = BestEvent(G, 0)
                        RelY = BestEvent(G, 1)
                        'replace these with routines that reconstruct the sequences in a way that take sinto account simmultaneous
                        'movement of separate pieces
                        Call ModSeqNum(0)
                        Call MakeTreeSeqs(ST, EN)
                        
                        P1 = Seq1
                        P2 = Seq2
                        d = Seq3
                        'Work out the maximum identity of parental sequences at the time when the recombination event occured
                        
                        
                        Call GetAge(G, Daught(), G, d, P1, P2, FMat(), SMat(), PermDiffs(), PermValid, AgeEvent())
                        
                        Call UnModSeqNum(0)
                    End If
                End If
        
                
            End If
            
            Form1.SSPanel1.Caption = Trim(Str(G)) + " of " + Trim(Str(SEventNumber)) + " events mapped"
            Form1.SSPanel1.Refresh
            Form1.ProgressBar1 = (G / SEventNumber) * 70
            Form1.Refresh
        Next G
        
        
        
        Seq1 = xSeq1
        Seq2 = xSeq2
        Seq3 = xSeq3
        RelX = xrelx
        RelY = xrely
        NJFlag = xnjflag
    
    End If


Scan2Flag = 1
Do While CE <= SEventNumber
    ReDim GrpMaskSeq(Nextno + 100)
    oRelX = RelX
    oRelY = RelY
    RelX = 0: RelY = 0
    Call MakeBestEvent
    XOverList(0, 0).Eventnumber = XOverList(BestEvent(CE, 0), BestEvent(CE, 1)).Eventnumber
    
    Call ModSeqNum(0)
    XOverList(0, 0).Eventnumber = 0
    RelX = oRelX
    RelY = oRelY
    
    BE = XOverList(BestEvent(CE, 0), BestEvent(CE, 1)).Beginning '14711,1553
    EN = XOverList(BestEvent(CE, 0), BestEvent(CE, 1)).Ending
    'XX = SuperEventlist(XOverList(BestEvent(CE, 0), BestEvent(CE, 1)).Eventnumber)
    For Y = 0 To Nextno
        If Daught(CE, S2TraceBack(Y)) > 0 Then
            GrpMaskSeq(Y) = 1
        End If
    Next Y
    For z = LastNextNo + 1 To Nextno
        For Y = 1 To SEventNumber
            Daught(Y, z) = Daught(Y, S2TraceBack(z))
        Next Y
    Next z
    oSEventNumber = SEventNumber
    oEventNumber = Eventnumber
    Dim oSuperEventList() As Long
    ReDim oSuperEventList(Eventnumber)
    For X = 0 To Eventnumber
        oSuperEventList(X) = SuperEventList(X)
    Next X
    
    MaxDim = 0
    For X = 0 To PermNextNo
        If MaxDim < CurrentXOver(X) Then MaxDim = CurrentXOver(X)
    Next X
    ReDim BakXOList(PermNextNo, MaxDim), BakCurXOver(PermNextNo)
    For X = 0 To PermNextNo
    
       For Y = 1 To CurrentXOver(X)
            If SuperEventList(XOverList(X, Y).Eventnumber) <> CE Then
                BakCurXOver(X) = BakCurXOver(X) + 1
                BakXOList(X, BakCurXOver(X)) = XOverList(X, Y)
            End If
        Next Y
        
    Next X
    MaxDim = 0
    For X = 0 To PermNextNo
        If MaxDim < BCurrentXoverMi(X) Then MaxDim = BCurrentXoverMi(X)
    Next X
    ReDim BakXOlistMi(PermNextNo, MaxDim), BakCurXOverMi(PermNextNo)
    For X = 0 To PermNextNo
        
        For Y = 1 To BCurrentXoverMi(X)
            If SuperEventList(BestXOListMi(X, Y).Eventnumber) <> CE Then
                BakCurXOverMi(X) = BakCurXOverMi(X) + 1
                BakXOlistMi(X, BakCurXOverMi(X)) = BestXOListMi(X, Y)
            End If
        Next Y
    Next X
    MaxDim = 0
    For X = 0 To PermNextNo
        If MaxDim < BCurrentXoverMa(X) Then MaxDim = BCurrentXoverMa(X)
    Next X
    ReDim BakXOlistMa(PermNextNo, MaxDim), BakCurXOverMa(PermNextNo)
    For X = 0 To PermNextNo
        For Y = 1 To BCurrentXoverMa(X)
            If SuperEventList(BestXOListMa(X, Y).Eventnumber) <> CE Then
                BakCurXOverMa(X) = BakCurXOverMa(X) + 1
                BakXOlistMa(X, BakCurXOverMa(X)) = BestXOListMa(X, Y)
            End If
        Next Y
    Next X
    'daught, minpar,majpar,nopini, and dscores must be sorted out - ie take new scores and add to/replace old scores
    
    StepNo = 0
    Dim oBestEvent() As Long, oDaught() As Byte, oMinorPar() As Byte, oMajorPar() As Byte, oNOPINI() As Byte, oDScores() As Double
    'note that the orientation of the eliments is reversed in the odaugt, o minorpar etc.
    ReDim oDaught(Nextno, SEventNumber), oMinorPar(Nextno, SEventNumber), oMajorPar(Nextno, SEventNumber)
    ReDim oNOPINI(2, SEventNumber)
    ReDim oDScores(20, 2, SEventNumber)
    ReDim oBestEvent(1, SEventNumber)
    
    For Y = 1 To SEventNumber
        For z = 0 To Nextno
            oDaught(z, Y) = Daught(Y, S2TraceBack(z))
            oMinorPar(z, Y) = MinorPar(Y, S2TraceBack(z))
            oMajorPar(z, Y) = MajorPar(Y, S2TraceBack(z))
        Next z
        For z = 0 To 2
            oNOPINI(z, Y) = NOPINI(z, Y)
            For A = 0 To 20
                oDScores(A, z, Y) = DScores(A, z, Y)
            Next A
        Next z
        For z = 0 To 1
            oBestEvent(z, Y) = BestEvent(Y, z)
        Next z
    Next Y
    XOverListSize = 100
    ReDim CurrentXOver(Nextno), XOverList(Nextno, 100)
    ReDim SuperEventList(Nextno)
    
    Call DoRDP(0, 0)
    
    'Make a tempxolist with only the best event for each of the sequences is odaught(ce)
    'Make comparrison array for the current region
    Dim FirstCycle As Long, s(2) As Long, ReplaceE As Long, DoneEventX() As Byte, CompA() As Long, RSize(5) As Long, BPos2 As Long, EPos2 As Long, TMatch(1) As Double, BestMatch As Double, WinMatch(1) As Long
    
    
    
    ReDim CompA(Len(StrainSeq(0)))
    ReDim DoneEventX(SEventNumber)
    Dummy = MakeOLSeq(Len(StrainSeq(0)), BE, EN, RSize(0), CompA(0)) '4598,5975
    FirstCycle = 0
    XX = UBound(Daught, 1)
    
    Do
        BestMatch = 0: WinMatch(0) = 0: WinMatch(1) = 0
        For z = 0 To Nextno
            If oDaught(S2TraceBack(z), CE) > 0 Then
                For X = 0 To Nextno
                    For Y = 1 To CurrentXOver(X)
                        If DoneEventX(SuperEventList(XOverList(X, Y).Eventnumber)) = 0 Then
                            s(0) = S2TraceBack(XOverList(X, Y).Daughter)
                            s(1) = S2TraceBack(XOverList(X, Y).MajorP)
                            s(2) = S2TraceBack(XOverList(X, Y).MinorP)
                            For A = 0 To 2
                                If s(0) = S2TraceBack(z) Then
                                    Exit For
                                End If
                            Next A
                            If A < 3 Then
                                BPos2 = XOverList(X, Y).Beginning
                                EPos2 = XOverList(X, Y).Ending
                                OLSize = FindOverlap(Len(StrainSeq(0)), BPos2, EPos2, RSize(0), CompA(0))
                                If OLSize > 0 Then
                                    TMatch(1) = (OLSize * 2) / (RSize(0) + RSize(1))
                                Else
                                    TMatch(1) = 0
                                End If
                                           
                                TMatch(1) = CLng(TMatch(1) * 100000) / 100000
                                If TMatch(1) >= BestMatch Then
                                    BestMatch = TMatch(1)
                                    WinMatch(0) = X
                                    WinMatch(1) = Y
                                End If
                            End If
                        End If
                    Next Y
                Next X
            End If
        Next z
        'replace ce with supereventlist(xoverlist(winmatch(0),winmatch(1).eventnumber)
        'BakXOlist(PermNextNo, MaxDim), BakCurXOver(PermNextNo)
        If BestMatch = 0 Then Exit Do
        Dim WinE As Long
        WinE = SuperEventList(XOverList(WinMatch(0), WinMatch(1)).Eventnumber)
        DoneEventX(WinE) = 1
        If FirstCycle = 0 Then
            ReplaceE = CE
        ElseIf FirstCycle = 1 Then
            ReplaceE = oSEventNumber + 1
        Else
            ReplaceE = ReplaceE + 1
        End If
        If CE = 3 Then
            X = X
        End If
        
        Call AddToBak(SuperEventList(), ReplaceE, WinE, oEventNumber, S2TraceBack(), oSuperEventList(), BakXOlistMi(), BakCurXOverMi(), BestXOListMi(), BCurrentXoverMi())
        Call AddToBak(SuperEventList(), ReplaceE, WinE, oEventNumber, S2TraceBack(), oSuperEventList(), BakXOlistMa(), BakCurXOverMa(), BestXOListMa(), BCurrentXoverMa())
        Call AddToBak(SuperEventList(), ReplaceE, WinE, oEventNumber, S2TraceBack(), oSuperEventList(), BakXOList(), BakCurXOver(), XOverList(), CurrentXOver())

        If ReplaceE > UBound(oDaught, 2) Then
            ReDim Preserve oDaught(Nextno, ReplaceE + 100), oMinorPar(Nextno, ReplaceE + 100), oMajorPar(Nextno, ReplaceE + 100)
            ReDim Preserve oNOPINI(2, ReplaceE + 100)
            ReDim Preserve oDScores(20, 2, ReplaceE + 100)
            ReDim Preserve oBestEvent(1, ReplaceE + 100)
        End If
        For z = 0 To Nextno
            oDaught(z, ReplaceE) = Daught(WinE, z)
            oMinorPar(z, ReplaceE) = MinorPar(WinE, z)
            oMajorPar(z, ReplaceE) = MajorPar(WinE, z)
        Next z
        For z = 0 To 2
            oNOPINI(z, ReplaceE) = NOPINI(z, WinE)
            For A = 0 To 20
                oDScores(A, z, ReplaceE) = DScores(A, z, WinE)
            Next A
        Next z
        
        'oBestEvent(0, ReplaceE) = BestEvent(WinE, 0): oBestEvent(1, ReplaceE) = BestEvent(WinE, 1)
        If ReplaceE > UBound(AgeEvent, 2) Then
            ReDim Preserve AgeEvent(1, ReplaceE + 10)
        End If
        AgeEvent(0, ReplaceE) = AgeEvent(0, CE): AgeEvent(1, ReplaceE) = AgeEvent(1, CE)
        FirstCycle = FirstCycle + 1
        Exit Do
    Loop
    'BestXOListMiDim MaxDim As Long, BakXOlist() As XOverDefine, BakCurXOver() As Integer, BakXOlistMi() As XOverDefine, BakCurXOverMi() As Integer, BakXOlistMa() As XOverDefine, BakCurXOverMa() As Integer
    
    
    'got to erase event ce from bakxolists
    If FirstCycle = 0 Then
        Call EraseEvidence(Nextno, CE, BakXOList(), BakCurXOver(), oSuperEventList())
        Call EraseEvidence(Nextno, CE, BakXOlistMi(), BakCurXOverMi(), oSuperEventList())
        Call EraseEvidence(Nextno, CE, BakXOlistMa(), BakCurXOverMa(), oSuperEventList())
        
    End If
    
    
    'oEventNumber = Eventnumber + oEventNumber
    
    ReDim Preserve SuperEventList(UBound(oSuperEventList, 1))
    For X = 1 To UBound(oSuperEventList, 1)
        SuperEventList(X) = oSuperEventList(X)
    Next X
    
    For X = 1 To UBound(oSuperEventList, 1)
        If SuperEventList(X) > 0 Then Y = X
    Next X
    ReDim Preserve SuperEventList(Y)
    Eventnumber = Y
    If FirstCycle > 0 Then
        SEventNumber = oSEventNumber + FirstCycle - 1
    Else
        SEventNumber = oSEventNumber
    End If
    ReDim Daught(SEventNumber, Nextno), MinorPar(SEventNumber, Nextno), MajorPar(SEventNumber, Nextno)
    ReDim NOPINI(2, SEventNumber)
    ReDim DScores(20, 2, SEventNumber)
    ReDim BestEvent(SEventNumber, 1)
    For Y = 1 To SEventNumber
        For z = 0 To Nextno
            Daught(Y, S2TraceBack(z)) = oDaught(z, Y)
            MinorPar(Y, S2TraceBack(z)) = oMinorPar(z, Y)
            MajorPar(Y, S2TraceBack(z)) = oMajorPar(z, Y)
        Next z
        For z = 0 To 2
            NOPINI(z, Y) = oNOPINI(z, Y)
            For A = 0 To 20
                DScores(A, z, Y) = oDScores(A, z, Y)
            Next A
        Next z
        For z = 0 To 1
            BestEvent(Y, z) = oBestEvent(z, Y)
        Next z
    Next Y
    MaxDim = 0
    For X = 0 To PermNextNo
        If MaxDim < BakCurXOver(X) Then MaxDim = BakCurXOver(X)
    Next X
    ReDim XOverList(PermNextNo, MaxDim), CurrentXOver(PermNextNo)
    For X = 0 To PermNextNo
        CurrentXOver(X) = BakCurXOver(X)
        For Y = 1 To CurrentXOver(X)
                
                XOverList(X, Y) = BakXOList(X, Y)
                
        Next Y
        
    Next X
    MaxDim = 0
    For X = 0 To PermNextNo
        If MaxDim < BakCurXOverMi(X) Then MaxDim = BakCurXOverMi(X)
    Next X
    ReDim BestXOListMi(PermNextNo, MaxDim), BCurrentXoverMi(PermNextNo)
    For X = 0 To PermNextNo
        BCurrentXoverMi(X) = BakCurXOverMi(X)
        For Y = 1 To BCurrentXoverMi(X)
            
                BestXOListMi(X, Y) = BakXOlistMi(X, Y)
            
        Next Y
        
    Next X
    MaxDim = 0
    For X = 0 To PermNextNo
        If MaxDim < BakCurXOverMa(X) Then MaxDim = BakCurXOverMa(X)
    Next X
    ReDim BestXOListMa(PermNextNo, MaxDim), BCurrentXoverMa(PermNextNo)
    For X = 0 To PermNextNo
        BCurrentXoverMa(X) = BakCurXOverMa(X)
        For Y = 1 To BCurrentXoverMa(X)
            
                BestXOListMa(X, Y) = BakXOlistMa(X, Y)
            
        Next Y
        
    Next X
    'XX = XOverList(BestEvent(CE + 1, 0), BestEvent(CE + 1, 1)).Beginning
    'XX = XOverList(BestEvent(CE + 1, 0), BestEvent(CE + 1, 1)).Ending
    
    CE = CE + 1
    
    XX = UBound(BestEvent, 1)
Loop


Call MakeBestEvent



SelGrpFlag = oSelGrpFlag
For X = 0 To Nextno
    GrpMaskSeq(X) = oGrpMask(X)
Next X

End Sub
Public Sub DrawGMap()

Dim GMap() As Integer
ReDim GMap(Len(StrainSeq(0)), Nextno)

For X = 0 To Nextno
    For Y = 1 To Len(StrainSeq(0))
        If SeqNum(Y, X) = 46 Then
            'check for missing data
            Cnt = 0
            For z = Y To Y + 10 'len(strainseq(0))
                If z <= Len(StrainSeq(0)) Then
                    If SeqNum(z, X) = 46 Then
                        Cnt = Cnt + 1
                    End If
                Else
                    Exit For
                    Cnt = 10
                End If
            Next z
            If Cnt >= 10 Then
                For z = Y To Len(StrainSeq(0))
                    If SeqNum(z, X) = 46 Then
                        GMap(z, X) = -1
                    Else
                        Y = z
                        Exit For
                    End If
                Next z
            
            End If
        End If
        
    Next Y
Next X



For X = 0 To Nextno
    For Y = 1 To Len(StrainSeq(0))
        If GMap(Y, X) <> -1 Then
            If SeqNum(Y, X) = SeqNum(Y, 0) And SeqNum(Y, X) = SeqNum(Y, 1) Then
            
            ElseIf SeqNum(Y, X) = SeqNum(Y, 0) Then
                GMap(Y, X) = 1
            ElseIf SeqNum(Y, X) = SeqNum(Y, 1) Then
                GMap(Y, X) = 2
            Else
                GMap(Y, X) = 100
            End If
        End If
        
    Next Y
Next X

On Error Resume Next
Kill "tmp2.emf"
On Error GoTo 0
SEMFnameII = "tmp2.emf"
semfname$ = "tmp2.emf"

Dim SP As Long, EP As Long
Dim PColIn As Long, HFactor As Double, WFactor As Double
Dim OldFont As Long, OldPen As Long, Pen As Long, LOffset As Long, TOffset As Long, MhDC As Long, EMFCls As Long

Dim LPn As LOGPEN
Dim red As Long, green As Long, Green2 As Long, blue As Long
Dim rct As RECT
Dim LoFnt As Long

rct.left = 0
rct.Top = 0
rct.right = Len(StrainSeq(0)) * 20
rct.Bottom = 16000


HFactor = 200
WFactor = 0.5 '500 / Len(StrainSeq(0))
LOffset = 70
TOffset = 40
    
Form1.Picture1.AutoRedraw = False
MhDC = CreateEnhMetaFile(Form1.Picture1.hDC, semfname$, rct, "")
Form1.Picture1.AutoRedraw = True

'Get original Metafile font and pen
LoFnt = CreateFont(20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "Arial")
OldFont = SelectObject(MhDC, LoFnt)
Pen = CreatePenIndirect(LPn)
OldPen = SelectObject(MhDC, Pen)


    
    
'Draw lines
Dim PntAPI As POINTAPI
For X = 0 To Nextno
    For Y = 0 To Len(StrainSeq(0))
        If GMap(Y, X) <> 0 Then
            If GMap(Y, X) = -1 Then
                LPn.lopnColor = RGB(196, 196, 196)
                Pen = CreatePenIndirect(LPn)
                SelectObject MhDC, Pen
                dl = 3
            ElseIf GMap(Y, X) = 1 Then
                LPn.lopnColor = RGB(255, 0, 0)
                Pen = CreatePenIndirect(LPn)
                SelectObject MhDC, Pen
                dl = 5
            ElseIf GMap(Y, X) = 2 Then
                LPn.lopnColor = RGB(0, 0, 255)
                Pen = CreatePenIndirect(LPn)
                SelectObject MhDC, Pen
                dl = 5
            Else
                LPn.lopnColor = RGB(0, 0, 0)
                Pen = CreatePenIndirect(LPn)
                SelectObject MhDC, Pen
                dl = 10
            End If
            
            MoveToEx MhDC, LOffset + WFactor * Y, TOffset + TOffset * X - dl, PntAPI
            LineTo MhDC, LOffset + WFactor * Y, TOffset + TOffset * X + dl
        End If
    Next Y
Next X

LPn.lopnColor = RGB(0, 0, 0)
Pen = CreatePenIndirect(LPn)
SelectObject MhDC, Pen
'Draw backbones of sequences

For X = 0 To Nextno

    MoveToEx MhDC, LOffset, TOffset + TOffset * X, PntAPI
    LineTo MhDC, LOffset + WFactor * Len(StrainSeq(0)), TOffset + TOffset * X
    TextOut MhDC, LOffset + WFactor * Len(StrainSeq(0)) + 20, TOffset + TOffset * X - 10, StraiName(X), Len(StraiName(X))
Next X
 


'Clear up and close  emf
Pen = SelectObject(MhDC, OldPen)
DeleteObject (Pen)
LoFnt = SelectObject(MhDC, OldFont)
DeleteObject (LoFnt)
EMFCls = CloseEnhMetaFile(MhDC)
Dummy = DeleteEnhMetaFile(EMFCls)
    
Clipboard.Clear
Clipboard.SetData LoadPicture("tmp2.emf"), 3
On Error Resume Next
On Error GoTo 0

End Sub
Public Sub AlphabetiseNames()

Dim SCx() As Long
ReDim SCx(Nextno)
For X = 0 To Nextno
    
    If InStr(1, StraiName(X), "{", vbBinaryCompare) > 0 Or Len(StraiName(X)) < 5 Then
        
    Else
        
        Do While Val(left(StraiName(X), 1)) > 0 Or left(StraiName(X), 1) = "0"
            StraiName(X) = Mid(StraiName(X), 2, Len(StraiName(X)) - 1)
            If Mid(StraiName(X), 3, 1) = "_" Then Exit Do 'Or Mid(StraiName(X), 3, 1) = "."
        Loop
        
        
        X = X
    End If
    SCx(X) = 0
    SCx(X) = SCx(X) + Asc(Mid(StraiName(X), 1, 1))
    SCx(X) = SCx(X) * 1000
    SCx(X) = SCx(X) + Asc(Mid(StraiName(X), 2, 1)) * 100
    SCx(X) = SCx(X) + Asc(Mid(StraiName(X), 3, 1)) * 10
    If Len(StraiName(X)) > 3 Then
        SCx(X) = SCx(X) + Asc(Mid(StraiName(X), 4, 1))
    End If
Next X
WinX = 0
For X = 0 To Nextno
    MS = 0
    For Y = 0 To Nextno
        If SCx(Y) > MS Then
            MS = SCx(Y)
            WinX = Y
        End If
    Next Y
    If MS = 0 Then
        Exit For
    Else
        Print #1, ">" + StraiName(WinX)
        Print #1, StrainSeq(WinX)
        SCx(WinX) = 0
    End If
Next X

'Close #2
Close #1

End Sub
Public Sub EnableFrame6()
'Form3.Frame6(0).Visible = True
'Form3.Frame6(2).Visible = False
'Exit Sub
With Form3
    .Frame6(0).Enabled = True
    .Label1(45).Enabled = True
    .Label1(46).Enabled = True
    .Text1(33).Enabled = True
    .Text1(34).Enabled = True
    .Text1(33).BackColor = RGB(255, 255, 255)
    .Text1(34).BackColor = RGB(255, 255, 255)
    .Text1(33).ForeColor = 0
    .Text1(34).ForeColor = 0
End With


End Sub
Public Sub DisableFrame6()
'Form3.Frame6(0).Visible = False
'Form3.Frame6(2).Visible = True
'Exit Sub
With Form3
    .Frame6(0).Enabled = False
    .Label1(45).Enabled = False
    .Label1(46).Enabled = False
    .Text1(33).Enabled = False
    .Text1(34).Enabled = False
    .Text1(33).ForeColor = RGB(128, 128, 128)
    .Text1(34).ForeColor = RGB(128, 128, 128)
    .Text1(33).BackColor = Form1.BackColor
    .Text1(34).BackColor = Form1.BackColor
End With


End Sub
Public Sub SetMLModel()


If ModelTestFlag = 1 Then
    Form3.Label1(69) = "Automatic model selection"
    Form3.Label1(24).Enabled = False
    Form3.Command28(43).Enabled = False
    Form3.Command28(44).Enabled = False
    Form3.Label1(59).Enabled = False
    Form3.Label1(60).Enabled = False
    Form3.Label21(46).Enabled = False
    Form3.Label21(42).Enabled = False
    Form3.Label21(52).Enabled = False
    Form3.Text23(34).BackColor = Form1.BackColor
    Form3.Text23(42).BackColor = Form1.BackColor
    Form3.Text1(18).BackColor = Form1.BackColor
    Form3.Text1(17).BackColor = Form1.BackColor
    Form3.Text23(34).Enabled = False
    Form3.Text23(42).Enabled = False
    Form3.Text1(18).Enabled = False
    Form3.Text1(17).Enabled = False
    
        
Else
    Form3.Label1(69) = "User specified model"
    Form3.Label1(24).Enabled = True
    Form3.Command28(43).Enabled = True
    Form3.Command28(44).Enabled = True
    Form3.Label1(59).Enabled = True
    Form3.Label1(60).Enabled = True
    Form3.Label21(46).Enabled = True
    Form3.Label21(42).Enabled = True
    Form3.Label21(52).Enabled = True
    Form3.Text23(34).BackColor = QBColor(15)
    Form3.Text23(42).BackColor = QBColor(15)
    Form3.Text1(18).BackColor = QBColor(15)
    Form3.Text1(17).BackColor = QBColor(15)
    Form3.Text23(34).Enabled = True
    Form3.Text23(42).Enabled = True
    Form3.Text1(18).Enabled = True
    Form3.Text1(17).Enabled = True
    If TPModel = 0 Then
        Form3.Label1(24).Caption = "Jukes and Cantor, 1969"
        Form3.Frame21(6).Enabled = False
        
        'Disable Transition transversion ratio option
        Form3.Text1(17).Text = "0.5"
        Form3.Text1(17).Enabled = False
        Form3.Text1(17).BackColor = Form1.BackColor
        Form3.Label1(59).Enabled = False
        
        'Disable base frequency estimate option
        Form3.Label21(46).Enabled = False
        Form3.Command28(44).Enabled = False
        
        
        
        Form3.Text23(30).BackColor = Form1.BackColor
        Form3.Text23(31).BackColor = Form1.BackColor
        Form3.Text23(32).BackColor = Form1.BackColor
        Form3.Text23(33).BackColor = Form1.BackColor
        Form3.Text23(30).Enabled = False
        Form3.Text23(31).Enabled = False
        Form3.Text23(32).Enabled = False
        Form3.Text23(33).Enabled = False
        Form3.Label21(36).Enabled = False
        Form3.Label21(37).Enabled = False
        Form3.Label21(38).Enabled = False
        Form3.Label21(39).Enabled = False
    ElseIf TPModel = 1 Then
        Form3.Label1(24).Caption = "Kimura, 1980"
        
        'Enable transition transversion rate selection option
        Form3.Text1(17).Text = TPTVRat
        Form3.Text1(17).Enabled = True
        Form3.Text1(17).BackColor = QBColor(15)
        Form3.Label1(59).Enabled = True
        
        'Disable base frequency estimate option
        Form3.Label21(46).Enabled = False
        Form3.Command28(44).Enabled = False
        
        
        
    ElseIf TPModel = 2 Then
        Form3.Label1(24).Caption = "Felsenstein, 1981"
        'Disable Transition transversion ratio option
        Form3.Text1(17).Text = "0.5"
        Form3.Text1(17).Enabled = False
        Form3.Text1(17).BackColor = Form1.BackColor
        Form3.Label1(59).Enabled = False
        
        'Enable base frequency estimate option
        Form3.Label21(46).Enabled = True
        Form3.Command28(44).Enabled = True
        
        
        
    ElseIf TPModel = 3 Then
        Form3.Label1(24).Caption = "Felsenstein, 1984"
        'Enable transition transversion rate selection option
        Form3.Text1(17).Text = TPTVRat
        Form3.Text1(17).Enabled = True
        Form3.Text1(17).BackColor = QBColor(15)
        Form3.Label1(59).Enabled = True
        
        'Enable base frequency estimate option
        Form3.Label21(46).Enabled = True
        Form3.Command28(44).Enabled = True
        
        
        
       
        
    ElseIf TPModel = 4 Then
        Form3.Label1(24).Caption = "Tamura and Nei, 1993"
         'Enable transition transversion rate selection option
        Form3.Text1(17).Text = TPTVRat
        Form3.Text1(17).Enabled = True
        Form3.Text1(17).BackColor = QBColor(15)
        Form3.Label1(59).Enabled = True
        
        'Enable base frequency estimate option
        Form3.Label21(46).Enabled = True
        Form3.Command28(44).Enabled = True
        
        
        
    
    ElseIf TPModel = 6 Then
        Form3.Label1(24).Caption = "Hasagawa, Kishino and Yano, 1985"
         'Enable transition transversion rate selection option
        Form3.Text1(17).Text = TPTVRat
        Form3.Text1(17).Enabled = True
        Form3.Text1(17).BackColor = QBColor(15)
        Form3.Label1(59).Enabled = True
        
        'Enable base frequency estimate option
        Form3.Label21(46).Enabled = True
        Form3.Command28(44).Enabled = True
        
        
        Form3.Frame21(6).Enabled = True
        
    
        
    ElseIf TPModel = 5 Then
        Form3.Label1(24).Caption = "General Time Reversable"
        'Disable Transition transversion ratio option
        Form3.Text1(17).Text = "0.5"
        Form3.Text1(17).Enabled = False
        Form3.Text1(17).BackColor = Form1.BackColor
        Form3.Label1(59).Enabled = False
        
        'Enable base frequency estimate option
        Form3.Label21(46).Enabled = True
        Form3.Command28(44).Enabled = True
        
      
        
    End If
End If
End Sub

Public Sub SetTBGamma()
If TBGamma = 0 Then
    Form3.Label21(51).Caption = "No rate variation accross sites"
    Form3.Text1(20).Enabled = False
    Form3.Text1(20).ForeColor = RGB(128, 128, 128)
    Form3.Text1(20).BackColor = Form1.BackColor
    Form3.Label1(62).Enabled = False
ElseIf TBGamma = 1 Then

    Form3.Label21(51).Caption = "Gamma-distributed rate variation"
    
    Form3.Text1(20).Enabled = True
    Form3.Text1(20).ForeColor = 0
    Form3.Text1(20).BackColor = RGB(255, 255, 255)
    Form3.Label1(62).Enabled = True
    
ElseIf TBGamma = 2 Then
    Form3.Label21(51).Caption = "Auto-correlated gamma-distributed rate variation"
ElseIf TBGamma = 3 Then
ElseIf TBGamma = 4 Then
Else
End If
End Sub
Public Sub SetTBModel()

If TBModel = 0 Then
    Form3.Label1(61).Caption = "All 6 substitution types are equally likely"
ElseIf TBModel = 1 Then
    Form3.Label1(61).Caption = "Transitions/transversions can be unequally likely "
Else
    Form3.Label1(61).Caption = "All 6 substitution types can be unequally likely"
End If
 Form3.Label1(61).Refresh
 End Sub
Public Sub SetF3Vals(BMFlag)


With Form3
        OLSeq = .List1.TopIndex

        If SpacerFlag = 0 Then
            .Option2.Value = 1
        ElseIf SpacerFlag = 1 Then
            .Option3.Value = 1
        ElseIf SpacerFlag = 2 Then
            .Option4.Value = 1
        ElseIf SpacerFlag = 3 Then
            .Option5.Value = 1
        ElseIf SpacerFlag = 4 Then
            .Option6.Value = 1
        ElseIf SpacerFlag = 5 Then
            .Option6.Value = 1
        End If

        .Frame4(SpacerFlag).Visible = True
        .Frame4(SpacerFlag).ZOrder
        HomologyIndicatorT = HomologyIndicator
        CircularFlagT = CircularFlag
        WeightedFlagT = WeightedFlag
        ShowPlotFlagT = ShowPlotFlag
        MCFlagT = MCFlag
        MaDistanceX = MaDistance


        'do recrate settings
        xGCFlag = GCFlag
        .Text6(0) = StartRho
        .Text6(1) = BlockPen
        .Text6(2) = FreqCo
        .Text6(3) = FreqCoMD
        .Text6(4) = GCTractLen
        .Text6(5) = MCMCUpdates
        If BMFlag = 1 Then
        
            xBlockPen = BlockPen
            xStartRho = StartRho
            xMCMCUpdates = MCMCUpdates
            xFreqCo = FreqCo
            xFreqCoMD = FreqCoMD
            xGCFlag = GCFlag
            xGCTractLen = GCTractLen
            
            If GCFlag = 0 Then
                
                .Label2(4).Caption = "Do not use gene conversion model"
                .Label2(5).Enabled = False
                .Text6(4).Enabled = False
                .Text6(4).BackColor = Form1.BackColor
            Else
                .Label2(4).Caption = "Use gene conversion model"
                .Label2(5).Enabled = True
                .Text6(4).Enabled = True
                .Text6(4).BackColor = RGB(255, 255, 255)
            End If
            
            If CircularFlag = 1 Then
                .Label16 = "Sequences are circular"
            Else
                .Label16 = "Sequences are linear"
            End If
    
            .Text2.Text = Form1.Text5.Text
            .Text3.Text = Form1.Text1.Text
    
            If ShowPlotFlagT = 0 Then
                .Label20 = "Do not show plots during scan"
            ElseIf ShowPlotFlagT = 1 Then
                .Label20 = "Show plots during scan"
            ElseIf ShowPlotFlagT = 2 Then
                .Label20 = "Show overview during scan"
            End If
    
            If MCFlag = 1 Then
                .Label23 = "No multiple comparison correction"
            ElseIf MCFlag = 0 Then
                .Label23 = "Bonferroni correction"
            ElseIf MCFlag = 2 Then
                .Label23 = "Step down correction"
            End If
            xPermTypeFlag = PermTypeFlag
            .Text1(37) = GPerms
            If PermTypeFlag = 1 Then
            
                .Label1(54) = "Shuffle alignment columns"
            Else
                
                .Label1(54) = "Use SEQGEN parametric simulations"
            End If
            
               
    
            Dim TotT As Double
            xAllowConflict = AllowConflict
            If AllowConflict = 0 Then
                .Check13 = 1
            Else
                .Check13 = 0
            End If
            DontRefreshFlag = 1
           
            TotT = 0
            PNum = 0
            If DoScans(0, 0) = 1 Then
                .Check4.Value = 1
                TotT = TotT + AnalT(0)
                PNum = PNum + 1
            Else
                .Check4.Value = 0
            End If
     
            If DoScans(0, 1) = 1 Then
                .Check5.Value = 1
                TotT = TotT + AnalT(1)
                PNum = PNum + 1
            Else
                .Check5.Value = 0
            End If
    
            If DoScans(0, 2) = 1 Then
                .Check1.Value = 1
                TotT = TotT + AnalT(2)
                PNum = PNum + 1
            Else
                .Check1.Value = 0
            End If
    
            If DoScans(0, 3) = 1 Then
                .Check2.Value = 1
                TotT = TotT + AnalT(3)
                PNum = PNum + 1
            Else
                .Check2.Value = 0
            End If
    
            If DoScans(0, 4) = 1 Then
                .Check3.Value = 1
                TotT = TotT + AnalT(4)
                PNum = PNum + 1
            Else
                .Check3.Value = 0
            End If
    
            If DoScans(0, 5) = 1 Then
                .Check6.Value = 1
                TotT = TotT + AnalT(5)
                PNum = PNum + 1
            Else
                .Check6.Value = 0
            End If
            
            If DoScans(0, 8) = 1 Then
                .Check12.Value = 1
                TotT = TotT + AnalT(6)
                PNum = PNum + 1
            Else
                .Check12.Value = 0
            End If
                
            DontRefreshFlag = 0
            xConsensusProg = ConsensusProg
            
            
            If ConsensusProg = 0 Then
                .Label43 = "List all events"
            ElseIf ConsensusProg = 1 Then
                .Label43 = "List events detected by >1 method"
            ElseIf ConsensusProg = 2 Then
                .Label43 = "List events detected by >2 methods"
            ElseIf ConsensusProg = 3 Then
                .Label43 = "List events detected by >3 methods"
            ElseIf ConsensusProg = 4 Then
                .Label43 = "List events detected by >4 methods"
            ElseIf ConsensusProg = 5 Then
                .Label43 = "List events detected by >5 methods"
            ElseIf ConsensusProg = 6 Then
                .Label43 = "List events detected by >6 methods"
            End If
            
            .Check10 = ForcePhylE
            .Check9 = PolishBPFlag
            .Check11 = Realignflag
            
            Dim PWidth As Integer
            Dim CurXpos As Double
    
            PWidth = .Picture27.Width - 100
            CurXpos = 25
            .Picture27.Picture = LoadPicture()
   
            Call RefreshTimes
        End If
        .Picture27.Refresh
        
        
        
        If GCSeqTypeFlag = 0 Then
            .Label29.Caption = "Automatic detection of sequence type"
            .Label32.Enabled = False
            .Command28(4).Enabled = False
        ElseIf GCSeqTypeFlag = 1 Then
            .Label29.Caption = "Sequences are DNA"
            .Label32.Enabled = False
            .Command28(4).Enabled = False
        ElseIf GCSeqTypeFlag = 2 Then
            .Label29.Caption = "Sequences are DNA coding region"
            .Command28(4).Enabled = True
            .Label32.Enabled = True
        ElseIf GCSeqTypeFlag = 3 Then
            .Label29.Caption = "Sequences are protein"
            .Command28(4).Enabled = False
            .Label32.Enabled = False
        End If

        If GCIndelFlag = 0 Then
            .Label31.Caption = "Ignor indels"
        ElseIf GCIndelFlag = 1 Then
            .Label31.Caption = "Treat indel blocs as one polymorphism"
        ElseIf GCIndelFlag = 2 Then
            .Label31.Caption = "Treat each indel site as a polymorphism"
        End If

        If GCCodeFlag = 0 Then
            .Label32.Caption = "Use standard (nuclear) code"
        Else
            .Label32.Caption = "Use mammalian mitochondrial code"
        End If

        If GCMonoSiteFlag = 0 Then
            .Label39.Caption = "Do not use monomorphic sites"
        Else
            .Label39.Caption = "Use monomorphic sites"
        End If

        .Text21 = CStr(GCSeqRange(0))
        .Text22 = CStr(GCSeqRange(1))
        .Text9 = GCOutfileName

        If GCOutFlag = 0 Then
            .Label27.Caption = "Space separated output"
        ElseIf GCOutFlag = 1 Then
            .Label27.Caption = "Tab separated output"
        ElseIf GCOutFlag = 2 Then
            .Label27.Caption = "DIF-format spreadsheet output"
        ElseIf GCOutFlag = 3 Then
            .Label27.Caption = "Output in all formats"
        End If

        If GCOutFlagII = 0 Then
            .Label33.Caption = "Simple output"
            .Command29.Enabled = False
        Else
            .Label33.Caption = "Maximum output"
            'Command29.Enabled = True
        End If

        If GCSortFlag = 0 Then
            .Label34.Caption = "Sort fragment lists by P-Value"
        ElseIf GCSortFlag = 1 Then
            .Label34.Caption = "Sort lists alphabetically by name"
        ElseIf GCSortFlag = 2 Then
            .Label34.Caption = "Sort lists by P-Value then name"
        End If

        .Text10 = CStr(GCEndLen)
        .Text11 = CStr(GCOffsetAddjust)

        If GCLogFlag = 0 Then
            .Label28.Caption = "Write log file"
        ElseIf GCLogFlag = 1 Then
            .Label28.Caption = "Append existing log file"
        ElseIf GCLogFlag = 2 Then
            .Label28.Caption = "Do not write log file"
        End If

        .Text12 = CStr(GCMissmatchPen)
        .Text14 = CStr(GCMaxGlobFrags)
        .Text13 = CStr(GCMaxPairFrags)
        .Text15 = CStr(GCMinFragLen)
        .Text16 = CStr(GCMinPolyInFrag)
        .Text17 = CStr(GCMinPairScore)
        .Text18 = CStr(GCMaxOverlapFrags)
        .Text19 = CStr(GCNumPerms)
        .Text20 = GCMaxPermPVal

        If GCPermPolyFlag = 0 Then
            .Label59.Caption = "Use simple polymorphisms"
        Else
            .Label59.Caption = "Use only multiple polymorphisms"
        End If
        
        If GCtripletflag = 0 Then
            .Label39.Caption = "Scan sequence pairs"
            'Disbable large sections of the interface
            .Command28(2).Enabled = True
            .Command28(4).Enabled = True
            .Text21.Enabled = True
            .Text22.Enabled = True
            .Label48.Enabled = True
            .Label47.Enabled = True
            .Label32.Enabled = True
            .Label29.Enabled = True
            
            .Frame13.Enabled = True
            .Label38.Enabled = True
            .Text14.Enabled = True
            .Text13.Enabled = True
            .Label37.Enabled = True
            
            .Frame16.Enabled = True
            .Label26.Enabled = True
            .Label27.Enabled = True
            .Label33.Enabled = True
            .Label34.Enabled = True
            .Label30.Enabled = True
            .Label35.Enabled = True
            .Label28.Enabled = True
            .Label45.Enabled = True
            .Label46.Enabled = True
            .Label59.Enabled = True
            
            .Text19.Enabled = True
            .Text20.Enabled = True
            .Text9.Enabled = True
            .Text10.Enabled = True
            .Text11.Enabled = True
            .Command28(0).Enabled = True
            .Command28(5).Enabled = True
            .Command28(6).Enabled = True
            .Command28(1).Enabled = True
            .Command28(7).Enabled = True
            .Text21.BackColor = QBColor(15)
            .Text11.BackColor = QBColor(15)
            .Text10.BackColor = QBColor(15)
            .Text9.BackColor = QBColor(15)
            .Text20.BackColor = QBColor(15)
            .Text19.BackColor = QBColor(15)
            .Text14.BackColor = QBColor(15)
            .Text13.BackColor = QBColor(15)
            .Text22.BackColor = QBColor(15)
        Else
            .Label39.Caption = "Scan sequence triplets"
            'Disbable large sections of the interface
            .Command28(2).Enabled = False
            .Command28(4).Enabled = False
            .Text21.Enabled = False
            .Text22.Enabled = False
            .Label48.Enabled = False
            .Label47.Enabled = False
            .Label32.Enabled = False
            .Label29.Enabled = False
            
            .Frame13.Enabled = False
            .Label38.Enabled = False
            .Text14.Enabled = False
            .Text13.Enabled = False
            .Label37.Enabled = False
            
            .Frame16.Enabled = False
            .Label26.Enabled = False
            .Label27.Enabled = False
            .Label33.Enabled = False
            .Label34.Enabled = False
            .Label30.Enabled = False
            .Label35.Enabled = False
            .Label28.Enabled = False
            .Label45.Enabled = False
            .Label46.Enabled = False
            .Label59.Enabled = False
            
            .Text19.Enabled = False
            .Text20.Enabled = False
            .Text9.Enabled = False
            .Text10.Enabled = False
            .Text11.Enabled = False
            .Command28(0).Enabled = False
            .Command28(5).Enabled = False
            .Command28(6).Enabled = False
            .Command28(1).Enabled = False
            .Command28(7).Enabled = False
            .Text21.BackColor = Form1.Command1.BackColor
            .Text11.BackColor = Form1.Command1.BackColor
            .Text10.BackColor = Form1.Command1.BackColor
            .Text9.BackColor = Form1.Command1.BackColor
            .Text20.BackColor = Form1.Command1.BackColor
            .Text19.BackColor = Form1.Command1.BackColor
            .Text14.BackColor = Form1.Command1.BackColor
            .Text13.BackColor = Form1.Command1.BackColor
            .Text22.BackColor = Form1.Command1.BackColor
        End If
        
        xGCTripletFlag = GCtripletflag
        'store backed up GC variables
        xGCSeqTypeFlag = GCSeqTypeFlag
        xGCIndelFlag = GCIndelFlag
        xGCCodeFlag = GCCodeFlag
        xGCMonoSiteFlag = GCMonoSiteFlag
        ReDim xGCSeqRange(1)
        xGCSeqRange(0) = GCSeqRange(0)
        xGCSeqRange(1) = GCSeqRange(1)
        xGCOutfileName = GCOutfileName
        xGCOutFlag = GCOutFlag
        xGCOutFlagII = GCOutFlagII
        xGCSortFlag = GCSortFlag
        xGCEndLen = GCEndLen
        xGCOffsetAddjust = GCOffsetAddjust
        xGCLogFlag = GCLogFlag
        xGCMissmatchPen = GCMissmatchPen
        xGCMaxPairFrags = GCMaxPairFrags
        xGCMinFragLen = GCMinFragLen
        xGCMinPolyInFrag = GCMinPolyInFrag
        xGCMinPairScore = GCMinPairScore
        xGCMaxOverlapFrags = GCMaxOverlapFrags
        xGCNumPerms = GCNumPerms
        xGCMaxPermPVal = GCMaxPermPVal
        xGCPermPolyFlag = GCPermPolyFlag
        .Text1(0) = BSStepWin
        .Text1(1) = BSStepSize
        .Text1(5) = BSCutOff * 100
        .Text1(2) = BSBootReps
        .Text1(3) = BSRndNumSeed
        .Text1(28) = BSCoeffVar

        If BSSubModelFlag = 0 Then
            .Label1(4).Caption = "Jukes and Cantor, 1969"
            .Text1(4).Enabled = False
            .Text1(4).BackColor = Form1.BackColor
            .Text1(4).ForeColor = QBColor(8)
            .Text1(4) = 0.5
            .Label1(5).Enabled = False
            .Text1(28).Enabled = False
            .Text1(28).BackColor = Form1.BackColor
            .Text1(28).ForeColor = QBColor(8)
            .Label1(39).Enabled = False
            .Frame21(5).Enabled = False
            .Command28(24).Enabled = False

            For X = 26 To 29
                .Text23(X).BackColor = Form1.BackColor
                .Text23(X).Enabled = False
            Next 'X

            For X = 30 To 34
                .Label21(X).Enabled = False
            Next 'X
        ElseIf BSSubModelFlag = 4 Then
            .Label1(4).Caption = "Similarities"
            .Text1(4).Enabled = False
            .Text1(4).BackColor = Form1.BackColor
            .Text1(4).ForeColor = QBColor(8)
            .Text1(4) = 0.5
            .Label1(5).Enabled = False
            .Text1(28).Enabled = False
            .Text1(28).BackColor = Form1.BackColor
            .Text1(28).ForeColor = QBColor(8)
            .Label1(39).Enabled = False
            .Frame21(5).Enabled = False
            .Command28(24).Enabled = False

            For X = 26 To 29
                .Text23(X).BackColor = Form1.BackColor
                .Text23(X).Enabled = False
            Next 'X

            For X = 30 To 34
                .Label21(X).Enabled = False
            Next 'X
        ElseIf BSSubModelFlag = 1 Then
            .Label1(4).Caption = "Kimura, 1980"
            .Text1(4).Enabled = True
            .Text1(4).BackColor = QBColor(15)
            .Text1(4).ForeColor = QBColor(0)
            .Text1(4) = BSTTRatio
            .Label1(5).Enabled = True
            .Text1(28).Enabled = False
            .Text1(28).BackColor = Form1.BackColor
            .Text1(28).ForeColor = QBColor(8)
            .Label1(39).Enabled = False
            .Frame21(5).Enabled = False
            .Command28(24).Enabled = False

            For X = 26 To 29
                .Text23(X).BackColor = Form1.BackColor
                .Text23(X).Enabled = False
            Next 'X

            For X = 30 To 34
                .Label21(X).Enabled = False
            Next 'X

        ElseIf BSSubModelFlag = 2 Then
            .Label1(4).Caption = "Jin  and  Nei, 1990"
            .Text1(4).Enabled = True
            .Text1(4).BackColor = QBColor(15)
            .Text1(4).ForeColor = QBColor(0)
            .Text1(4) = BSTTRatio
            .Text1(28).Enabled = True
            .Text1(28).BackColor = QBColor(15)
            .Text1(28).ForeColor = QBColor(0)
            .Label1(39).Enabled = True
            .Label1(5).Enabled = True
            .Frame21(5).Enabled = False
            .Command28(24).Enabled = False

            For X = 26 To 29
                .Text23(X).BackColor = Form1.BackColor
                .Text23(X).Enabled = False
            Next 'X

            For X = 30 To 34
                .Label21(X).Enabled = False
            Next 'X

        ElseIf BSSubModelFlag = 3 Then
            .Label1(4).Caption = "Felsenstein, 1984"
            .Text1(4).Enabled = True
            .Text1(4).BackColor = QBColor(15)
            .Text1(4).ForeColor = QBColor(0)
            .Text1(4) = BSTTRatio
            .Text1(28).Enabled = False
            .Text1(28).BackColor = Form1.BackColor
            .Text1(28).ForeColor = QBColor(8)
            .Label1(39).Enabled = False
            .Label1(5).Enabled = True
            .Frame21(5).Enabled = True
            .Command28(24).Enabled = True
            .Label21(34).Enabled = True

            If BSFreqFlag = 0 Then

                For X = 26 To 29
                    .Text23(X).BackColor = Form1.BackColor
                    .Text23(X).Enabled = False
                Next 'X

                For X = 30 To 33
                    .Label21(X).Enabled = False
                Next 'X

            Else

                For X = 26 To 29
                    .Text23(X).BackColor = QBColor(15)
                    .Text23(X).Enabled = True
                Next 'X

                For X = 30 To 33
                    .Label21(X).Enabled = True
                Next 'X

            End If

        End If

        If BSFreqFlag = 0 Then
            .Label21(34).Caption = "Estimate from alignment"
        Else
            .Label21(34).Caption = "User defined"
        End If

        .Text23(28) = BSFreqA
        .Text23(27) = BSFreqC
        .Text23(26) = BSFreqG
        .Text23(29) = BSFreqT

        If BSTypeFlag = 0 Then
            .Label1(37) = "Use distances"
        ElseIf BSTypeFlag = 1 Then
            .Label1(37) = "Use UPGMAs"
        ElseIf BSTypeFlag = 2 Then
            .Label1(37) = "Use neighbour joining trees"
        ElseIf BSTypeFlag = 3 Then
            .Label1(37) = "Use least squares trees"
        ElseIf BSTypeFlag = 4 Then
            .Label1(37) = "Use maximum likelihood trees"
        End If

        
        
        If BSPValFlag = 0 Then
            .Label1(48) = "Use bootstrap value as P-value"
        ElseIf BSPValFlag = 1 Then
            .Label1(48) = "Calculate binomial P-value"
        ElseIf BSPValFlag = 2 Then
            .Label1(48) = "Calculate Chi square P-value"
        End If
        
        .Text1(6) = BSStepWin
        .Text1(7) = BSStepSize
        .Text1(10) = BSBootReps
        '.Check7.Value = BSCCenterFlag
        .Check8.Value = BSCDecreaseStepFlag
        .Text1(8) = BSCDStepSize
        .Text1(11) = BSStepWin
        '.Text1(12) = BSCDBootReps
        .Text1(9) = BSCDSpan

        If .Check8.Value = 1 Then
            .Text1(8).Enabled = True
            .Text1(9).Enabled = True
            .Text1(11).Enabled = True
            .Text1(12).Enabled = True
            .Frame18.Enabled = True
            .Text1(8).BackColor = QBColor(15)
            .Text1(9).BackColor = QBColor(15)
            .Text1(8).ForeColor = QBColor(0)
            .Text1(9).ForeColor = QBColor(0)
            .Label1(7).Enabled = True
            .Label1(10).Enabled = True
            .Label1(14).Enabled = True
            .Label1(15).Enabled = True
        Else
            .Frame18.Enabled = False
            .Text1(8).Enabled = False
            .Text1(9).Enabled = False
            .Text1(11).Enabled = False
            .Text1(12).Enabled = False
            .Text1(8).BackColor = Form1.BackColor
            .Text1(9).BackColor = Form1.BackColor
            .Text1(8).ForeColor = QBColor(8)
            .Text1(9).ForeColor = QBColor(8)
            .Text1(11).BackColor = Form1.BackColor
            .Text1(12).BackColor = Form1.BackColor
            .Text1(11).ForeColor = QBColor(8)
            .Text1(12).ForeColor = QBColor(8)
            .Label1(7).Enabled = False
            .Label1(14).Enabled = False
            .Label1(15).Enabled = False
            .Label1(10).Enabled = False
        End If

        xBSStepWin = BSStepWin
        xBSStepSize = BSStepSize
        xBSCutoff = BSCutOff
        xBSBootReps = BSBootReps
        xBSRndNumSeed = BSRndNumSeed
        xBSSubModelFlag = BSSubModelFlag
        xBSTTRatio = BSTTRatio
        xBSStepWin = BSStepWin
        xBSStepSize = BSStepSize
        xBSBootReps = BSBootReps
        xMatPermNo = MatPermNo
        xBSCDecreaseStepFlag = BSCDecreaseStepFlag
        xBSCDStepSize = BSCDStepSize
        xBSStepWin = BSStepWin
        xBSCDSpan = BSCDSpan
        xBSTypeFlag = BSTypeFlag
        xBSFreqFlag = BSFreqFlag
        xMCProportionFlag = MCProportionFlag
        xMCWinFract = MCWinFract
        xMCStripGapsFlag = MCStripGapsFlag
        xBSPValFlag = BSPValFlag
        'store MaxChi variables

        If MCProportionFlag = 0 Then
            .Label1(11) = "Set window size"
            .Label1(27) = "# Variable sites per window"
            .Text1(21) = MCWinSize
        Else
            .Label1(11) = "Variable window size"
            .Label1(27) = "Fraction of variable sites per window"
            'If Nextno > 0 Then
            .Text1(21) = MCWinFract
            'End If
        End If

        If MCTripletFlag = 0 Then
            .Label1(13) = "Scan sequence triplets"
            '.Command28(14).Enabled = True
            '.Label1(16).ForeColor = QBColor(0)
        Else
            .Label1(13) = "Scan sequence pairs"
            '.Command28(14).Enabled = False
            '.Label1(16).Caption = "Use gaps"
            '.Label1(16).ForeColor = QBColor(8)
        End If

        If MCStripGapsFlag = 0 Then
            .Label1(16) = "Use gaps"
        Else
            .Label1(16) = "Strip gaps"
        End If
        
        
        
        
        
                
        
        xMCTripletFlag = MCTripletFlag
        xMCWinSize = MCWinSize
        xMCSteplen = MCSteplen
        xMCStart = MCStart
        xMCEnd = MCEnd
       
        
        
        
        xCWinFract = CWinFract
        xCProportionFlag = CProportionFlag
        xCWinSize = CWinSize
        
        
        If CProportionFlag = 0 Then
            .Label1(55) = "Set window size"
            .Label1(56) = "# Variable sites per window"
            .Text1(36) = CWinSize
        Else
            .Label1(55) = "Variable window size"
            .Label1(56) = "Fraction of variable sites per window"
            'If Nextno > 0 Then
            .Text1(36) = CWinFract
            'End If
        End If
        
        
        
        
        xLRDModel = LRDModel
        
        
        xSSGapFlag = SSGapFlag
        xSSVarPFlag = SSVarPFlag
        xSSOutlyerFlag = SSOutlyerFlag
        xSSRndSeed = SSRndSeed
        xSSWinLen = SSWinLen
        xSSStep = SSStep
        xSSNumPerms = SSNumPerms
        xSSNumPerms2 = SSNumPerms2
        xSSFastFlag = SSFastFlag
        
        .Text24(0) = SSWinLen
        .Text24(1) = SSStep
        .Text24(2) = SSNumPerms
        .Text24(3) = SSRndSeed
        .Text24(4) = SSNumPerms2
                
        
        
        

        If SSGapFlag = 0 Then
            .Label22(3) = "Strip gaps"
        Else
            .Label22(3) = "Use gaps"
        End If
        
        If SSVarPFlag = 0 Then
            .Label22(4) = "Use all positions"
        ElseIf SSVarPFlag = 1 Then
            .Label22(4) = "Use only 1/2/3/4 variable positions"
        Else
            .Label22(4) = "Use only 1/2/3 variable positions"
        End If
        
        If SSOutlyerFlag = 0 Then
            .Label22(5) = "Use randomised sequence"
        ElseIf SSOutlyerFlag = 1 Then
            .Label22(5) = "Use nearest outlyer"
        Else
            .Label22(5) = "Use most divergent sequence "
        End If
        
        If SSFastFlag = 0 Then
            .Label22(6) = "Do slow exhaustive scan"
        Else
            .Label22(6) = "Do fast scan"
        End If
        
        
        'Phylpro Options
        'xPPWinLen = PPWinLen
        xPPStripGaps = PPStripGaps
        xIncSelf = IncSelf
        'xPPSeed = PPSeed
        'xPPPerms = PPPerms
        .Text4(0) = PPWinLen
        .Text4(1) = PPSeed
        .Text4(2) = PPPerms
        
        If IncSelf = 1 Then
            .Label14(2) = "Use self comparrisons"
        Else
            .Label14(2) = "Do not use self comparrisons"
        End If
        
        If PPStripGaps = 0 Then
            .Label14(1) = "Ignor gaps"
        ElseIf PPStripGaps = 1 Then
            .Label14(1) = "Strip gaps"
        ElseIf PPStripGaps = 2 Then
            .Label14(1) = "Use gaps as fith character"
        End If
        
        If LRDModel = 0 Or LRDModel = 1 Then
            .Frame21(2).Enabled = False

            For X = 14 To 19
                .Label21(X).Enabled = False
            Next 'X

            .Label21(12).Enabled = True
            .Text23(11).Enabled = True
            .Text23(11).ForeColor = QBColor(0)
            .Text23(11).BackColor = QBColor(15)

            For X = 12 To 16
                .Text23(X).ForeColor = QBColor(8)
                .Text23(X).BackColor = Form3.BackColor
            Next 'X

        End If

        If LRDModel = 0 Then
            .Label21(0).Caption = "Hasegawa, Kishino and Yano, 1985"
        ElseIf LRDModel = 1 Then
            .Label21(0).Caption = "Falsenstein, 1984"
        ElseIf LRDModel = 2 Then
            .Label21(0).Caption = "Reversible process"
            .Frame21(2).Enabled = True

            For X = 14 To 19
                .Label21(X).Enabled = True
            Next 'X

            .Label21(12).Enabled = False
            .Text23(11).Enabled = False
            .Text23(11).ForeColor = QBColor(8)
            .Text23(11).BackColor = Form3.BackColor

            For X = 12 To 16
                .Text23(X).ForeColor = QBColor(0)
                .Text23(X).BackColor = QBColor(15)
            Next 'X

        End If

        .Text23(5) = LRDCategs
        .Text23(6) = LRDShape
        .Text23(11) = LRDTvRat
        .Text23(13) = LRDACCoeff
        .Text23(14) = LRDAGCoeff
        .Text23(15) = LRDATCoeff
        .Text23(12) = LRDCGCoeff
        .Text23(16) = LRDCTCoeff
        .Text23(17) = "1"
        xLRDBaseFreqFlag = LRDBaseFreqFlag

        If LRDModel = 1 Then
            .Command28(11).Enabled = True
            .Label21(13).Enabled = True
        Else
            .Command28(11).Enabled = False
            .Label21(13).Enabled = False
        End If

        If LRDBaseFreqFlag = 0 And LRDModel = 1 Then
            .Label21(13).Caption = "Estimate from alignment"

            For X = 8 To 11
                .Label21(X).Enabled = False
            Next 'X

            For X = 7 To 10
                .Text23(X).ForeColor = QBColor(8)
                .Text23(X).BackColor = Form3.BackColor
                .Text23(X).Enabled = False
            Next 'X

        ElseIf LRDBaseFreqFlag = 1 Then
            .Label21(13).Caption = "User defined"

            For X = 8 To 11
                .Label21(X).Enabled = True
            Next 'X

            For X = 7 To 10
                .Text23(X).ForeColor = QBColor(0)
                .Text23(X).BackColor = QBColor(15)
            Next 'X

        End If

        .Text23(9) = LRDAFreq
        .Text23(8) = LRDCFreq
        .Text23(7) = LRDGFreq
        .Text23(10) = LRDTFreq
        .Text23(2) = LRDCodon1
        .Text23(3) = LRDCodon2
        .Text23(4) = LRDCodon3
        .Text23(0) = LRDStep
        .Text23(1) = LRDWinLen
        xLRDRegion = LRDRegion
        xLRDWin = LRDWin
        If LRDRegion > 2 Then LRDRegion = 2
        If LRDRegion < 1 Then LRDRegion = 1
        If LRDRegion = 1 Then
            .Label21(2) = "Test one breakpoint"
            If LRDWin = 0 Then
                .Label21(40) = "Moving partition scan"
            End If
        Else
            .Label21(2) = "Test two breakpoints"
            If LRDWin = 0 Then
                .Label21(40) = "Moving partitions scan"
            End If
        End If
        
        If LRDWin = 0 Then
            .Text23(1).Enabled = False
            .Label21(41).Enabled = False
            .Text23(1).BackColor = Form1.BackColor
            '.Text23(1).ForeColor = 0
            .Label21(2).Enabled = True
            .Command28(40).Enabled = True
        Else
            .Label21(40) = "Sliding window scan"
            .Text23(1).Enabled = True
            .Text23(1).BackColor = QBColor(15)
            '.Text23(1).ForeColor = 0
            .Label21(41).Enabled = True
            .Label21(2).Enabled = False
            .Command28(40).Enabled = False
        End If
        
        'Reticulate stuff

        If RetSiteFlag = 1 Then
            .Label1(17) = "Use only binary sites"
        ElseIf RetSiteFlag = 2 Then
            .Label1(17) = "Use as transition/transversions"
        ElseIf RetSiteFlag = 3 Then
            .Label1(17) = "Use all sites"
        End If
        .Text5(0) = MatPermNo
        .Text5(1) = MatWinSize
        If Nextno > 0 Then
            .Label1(22).Enabled = True
            .Combo2.Enabled = True
        Else
            .Label1(22).Enabled = False
            .Combo2.Enabled = False
        End If
        xMatPermNo = MatPermNo
        xMatWinSize = MatWinSize
        If Nextno > 0 Then
            .Combo2.ListIndex = TypeSeqNumber
        End If
        If .Combo1.ListIndex = 0 Then
            .Command28(15).Enabled = True
            .Label1(17).Enabled = True
            .Label1(22).Enabled = False
            .Combo2.Enabled = False
            .Text5(0).Enabled = True
            .Text5(1).Enabled = False
            .Label1(26).Enabled = True
            .Label1(57).Enabled = False
        ElseIf .Combo1.ListIndex = 1 Then
            .Command28(15).Enabled = False
            .Label1(17).Enabled = False
            .Label1(22).Enabled = True
            .Combo2.Enabled = True
            .Text5(0).Enabled = False
            .Text5(1).Enabled = True
            .Label1(26).Enabled = False
            .Label1(57).Enabled = True
        ElseIf .Combo1.ListIndex = 2 Then
            .Command28(15).Enabled = False
            .Label1(17).Enabled = False
            .Label1(22).Enabled = True
            .Combo2.Enabled = True
            .Text5(0).Enabled = False
            .Text5(1).Enabled = False
            .Label1(26).Enabled = False
            .Label1(57).Enabled = False
        ElseIf .Combo1.ListIndex = 3 Then
            .Command28(15).Enabled = False
            .Label1(17).Enabled = False
            .Label1(22).Enabled = True
            .Combo2.Enabled = True
            .Text5(0).Enabled = False
            .Text5(1).Enabled = False
            .Label1(26).Enabled = False
            .Label1(57).Enabled = False
        ElseIf .Combo1.ListIndex = 4 Then
            .Command28(15).Enabled = False
            .Label1(17).Enabled = False
            .Label1(22).Enabled = True
            .Combo2.Enabled = True
            .Text5(0).Enabled = False
            .Text5(1).Enabled = True
            .Label1(26).Enabled = False
            .Label1(57).Enabled = True
        ElseIf .Combo1.ListIndex = 4 Then
            .Command28(15).Enabled = False
            .Label1(17).Enabled = False
            .Label1(22).Enabled = True
            .Combo2.Enabled = True
            .Text5(0).Enabled = True
            .Text5(1).Enabled = True
            .Label1(26).Enabled = True
            .Label1(57).Enabled = True
        End If
        Call DoColourScale
        If Nextno = 0 Then
            .Label1(22).Enabled = False
            .Combo2.Enabled = False
        End If
        

        'distance plot stuff
        xDPBFreqFlag = DPBFreqFlag
        xDPModelFlag = DPModelFlag
        xDPTVRatio = DPTVRatio
        .Text1(13) = DPWindow
        .Text1(14) = DPStep

        If DPBFreqFlag = 0 Then
            .Label21(20).Caption = "Estimate from alignment"
        Else
            .Label21(20).Caption = "User defined"
        End If

        .Text23(19) = DPBFreqA
        .Text23(20) = DPBFreqC
        .Text23(21) = DPBFreqG
        .Text23(18) = DPBFreqT
        .Text1(29) = DPCoeffVar

        If DPModelFlag = 0 Then
            .Label1(21).Caption = "Jukes and Cantor, 1969"
            .Frame21(3).Enabled = False
            tdbl = 0.5
            .Text1(15).Text = tdble
            .Text1(15).Enabled = False
            .Text1(15).BackColor = Form1.BackColor
            .Label1(20).Enabled = False
            .Text1(29).Enabled = False
            .Text1(29).BackColor = Form1.BackColor
            .Text1(29).ForeColor = QBColor(8)
            .Label1(40).Enabled = False
            .Label21(20).Enabled = False
            .Command28(17).Enabled = False
            .Text23(18).BackColor = Form1.BackColor
            .Text23(19).BackColor = Form1.BackColor
            .Text23(20).BackColor = Form1.BackColor
            .Text23(21).BackColor = Form1.BackColor
            .Text23(18).Enabled = False
            .Text23(19).Enabled = False
            .Text23(20).Enabled = False
            .Text23(21).Enabled = False
            .Label21(21).Enabled = False
            .Label21(22).Enabled = False
            .Label21(23).Enabled = False
            .Label21(24).Enabled = False
        ElseIf DPModelFlag = 4 Then
            .Label1(21).Caption = "Similarities"
            .Frame21(3).Enabled = False
            tdbl = 0.5
            .Text1(15).Text = tdbl
            .Text1(15).Enabled = False
            .Text1(15).BackColor = Form1.BackColor
            .Label1(20).Enabled = False
            .Text1(29).Enabled = False
            .Text1(29).BackColor = Form1.BackColor
            .Text1(29).ForeColor = QBColor(8)
            .Label1(40).Enabled = False
            .Label21(20).Enabled = False
            .Command28(17).Enabled = False
            .Text23(18).BackColor = Form1.BackColor
            .Text23(19).BackColor = Form1.BackColor
            .Text23(20).BackColor = Form1.BackColor
            .Text23(21).BackColor = Form1.BackColor
            .Text23(18).Enabled = False
            .Text23(19).Enabled = False
            .Text23(20).Enabled = False
            .Text23(21).Enabled = False
            .Label21(21).Enabled = False
            .Label21(22).Enabled = False
            .Label21(23).Enabled = False
            .Label21(24).Enabled = False
        ElseIf DPModelFlag = 1 Then
            .Label1(21).Caption = "Kimura, 1980"
            .Text1(15).Text = DPTVRatio
            .Text1(15).Enabled = True
            .Text1(15).BackColor = QBColor(15)
            .Label1(20).Enabled = True
            .Text1(29).Enabled = False
            .Text1(29).BackColor = Form1.BackColor
            .Text1(29).ForeColor = QBColor(8)
            .Label1(40).Enabled = False
            .Frame21(3).Enabled = False
            .Label21(20).Enabled = False
            .Command28(17).Enabled = False
            .Text23(18).BackColor = Form1.BackColor
            .Text23(19).BackColor = Form1.BackColor
            .Text23(20).BackColor = Form1.BackColor
            .Text23(21).BackColor = Form1.BackColor
            .Text23(18).Enabled = False
            .Text23(19).Enabled = False
            .Text23(20).Enabled = False
            .Text23(21).Enabled = False
            .Label21(21).Enabled = False
            .Label21(22).Enabled = False
            .Label21(23).Enabled = False
            .Label21(24).Enabled = False
        ElseIf DPModelFlag = 2 Then
            .Label1(21).Caption = "Jin  and  Nei, 1990"
            .Text1(15).Text = DPTVRatio
            .Text1(15).Enabled = True
            .Text1(15).BackColor = QBColor(15)
            .Label1(20).Enabled = True
            .Text1(29).Enabled = True
            .Text1(29).BackColor = QBColor(15)
            .Text1(29).ForeColor = QBColor(0)
            .Label1(40).Enabled = True
            .Frame21(3).Enabled = False
            .Label21(20).Enabled = False
            .Command28(17).Enabled = False
            .Text23(18).BackColor = Form1.BackColor
            .Text23(19).BackColor = Form1.BackColor
            .Text23(20).BackColor = Form1.BackColor
            .Text23(21).BackColor = Form1.BackColor
            .Text23(18).Enabled = False
            .Text23(19).Enabled = False
            .Text23(20).Enabled = False
            .Text23(21).Enabled = False
            .Label21(21).Enabled = False
            .Label21(22).Enabled = False
            .Label21(23).Enabled = False
            .Label21(24).Enabled = False
        ElseIf DPModelFlag = 3 Then
            .Label1(21).Caption = "Felsenstein, 1984"
            .Text1(15).Text = DPTVRatio
            .Text1(15).Enabled = True
            .Text1(15).BackColor = QBColor(15)
            .Label1(20).Enabled = True
            .Text1(29).Enabled = False
            .Text1(29).BackColor = Form1.BackColor
            .Text1(29).ForeColor = QBColor(8)
            .Label1(40).Enabled = False
            .Frame21(3).Enabled = True
            .Label21(20).Enabled = True
            .Command28(17).Enabled = True
            .Text23(18).BackColor = QBColor(15)
            .Text23(19).BackColor = QBColor(15)
            .Text23(20).BackColor = QBColor(15)
            .Text23(21).BackColor = QBColor(15)

            If DPBFreqFlag = 1 Then
                .Text23(18).Enabled = True
                .Text23(19).Enabled = True
                .Text23(20).Enabled = True
                .Text23(21).Enabled = True
                .Label21(21).Enabled = True
                .Label21(22).Enabled = True
                .Label21(23).Enabled = True
                .Label21(24).Enabled = True
            Else
                .Text23(18).Enabled = False
                .Text23(19).Enabled = False
                .Text23(20).Enabled = False
                .Text23(21).Enabled = False
                .Label21(21).Enabled = False
                .Label21(22).Enabled = False
                .Label21(23).Enabled = False
                .Label21(24).Enabled = False
            End If

        End If

        'TOPAL stuff
        '    Public TONumSeqs As Integer, TOPFlag As Integer, TORndNum As Integer, TOModel As Integer, TOTreeType As Integer, TOPerms As Integer, TOWinLen As Integer, TOStepSize As Integer, TOSmooth As Integer
        'Public TOHigh As Double, MatAverage As Double, TOPower As Double, TOPValCOff As Double, TOTvTs As Double, TOFreqFlag As Double, TOFreqA As Double, TOFreqC As Double, TOFreqG As Double, TOFreqT As Double
        '          xdpfreqflag = DPBFreqFlag
        '          xDPModelFlag = DPModelFlag
        '          xDPTVRatio = DPTVRatio
        .Text1(22) = TOWinLen
        .Text1(19) = TOStepSize

        If TOFreqFlag = 0 Then
            .Label21(29).Caption = "Estimate from alignment"
        Else
            .Label21(29).Caption = "User defined"
        End If

        .Text23(24) = TOFreqA
        .Text23(23) = TOFreqC
        .Text23(22) = TOFreqG
        .Text23(25) = TOFreqT
        .Text1(30) = TOCoeffVar

        If TOModel = 0 Then
            .Label1(25).Caption = "Jukes and Cantor, 1969"
            .Frame21(4).Enabled = False
            .Text1(16).Text = "0.5"
            .Text1(16).Enabled = False
            .Text1(16).BackColor = Form1.BackColor
            .Label1(28).Enabled = False
            .Text1(30).Enabled = False
            .Text1(30).BackColor = Form1.BackColor
            .Text1(30).ForeColor = QBColor(8)
            .Label1(41).Enabled = False
            .Label21(29).Enabled = False
            .Command28(20).Enabled = False
            .Text23(25).BackColor = Form1.BackColor
            .Text23(24).BackColor = Form1.BackColor
            .Text23(23).BackColor = Form1.BackColor
            .Text23(22).BackColor = Form1.BackColor
            .Text23(25).Enabled = False
            .Text23(24).Enabled = False
            .Text23(23).Enabled = False
            .Text23(22).Enabled = False
            .Label21(28).Enabled = False
            .Label21(27).Enabled = False
            .Label21(26).Enabled = False
            .Label21(25).Enabled = False
        ElseIf TOModel = 1 Then
            .Label1(25).Caption = "Kimura, 1980"
            .Text1(16).Text = TOTvTs
            .Text1(16).Enabled = True
            .Text1(16).BackColor = QBColor(15)
            .Label1(29).Enabled = True
            .Text1(30).Enabled = False
            .Text1(30).BackColor = Form1.BackColor
            .Text1(30).ForeColor = QBColor(8)
            .Label1(41).Enabled = False
            .Frame21(4).Enabled = False
            .Label21(29).Enabled = False
            .Command28(20).Enabled = False
            .Text23(25).BackColor = Form1.BackColor
            .Text23(24).BackColor = Form1.BackColor
            .Text23(23).BackColor = Form1.BackColor
            .Text23(22).BackColor = Form1.BackColor
            .Text23(25).Enabled = False
            .Text23(24).Enabled = False
            .Text23(23).Enabled = False
            .Text23(22).Enabled = False
            .Label21(28).Enabled = False
            .Label21(27).Enabled = False
            .Label21(26).Enabled = False
            .Label21(25).Enabled = False
        ElseIf TOModel = 2 Then
            .Label1(25).Caption = "Jin  and  Nei, 1990"
            .Text1(16).Text = TOTvTs
            .Text1(16).Enabled = True
            .Text1(16).BackColor = QBColor(15)
            .Label1(29).Enabled = True
            .Text1(30).Enabled = True
            .Text1(30).BackColor = QBColor(15)
            .Text1(30).ForeColor = QBColor(0)
            .Label1(41).Enabled = True
            .Frame21(4).Enabled = False
            .Label21(29).Enabled = False
            .Command28(20).Enabled = False
            .Text23(25).BackColor = Form1.BackColor
            .Text23(24).BackColor = Form1.BackColor
            .Text23(23).BackColor = Form1.BackColor
            .Text23(22).BackColor = Form1.BackColor
            .Text23(25).Enabled = False
            .Text23(24).Enabled = False
            .Text23(23).Enabled = False
            .Text23(22).Enabled = False
            .Label21(28).Enabled = False
            .Label21(27).Enabled = False
            .Label21(26).Enabled = False
            .Label21(25).Enabled = False
        ElseIf TOModel = 3 Then
            .Label1(25).Caption = "Felsenstein, 1984"
            .Text1(16).Text = TOTvTs
            .Text1(16).Enabled = True
            .Text1(16).BackColor = QBColor(15)
            .Label1(29).Enabled = True
            .Text1(30).Enabled = False
            .Text1(30).BackColor = Form1.BackColor
            .Text1(30).ForeColor = QBColor(8)
            .Label1(41).Enabled = False
            .Frame21(4).Enabled = True
            .Label21(29).Enabled = True
            .Command28(20).Enabled = True
            .Text23(25).BackColor = QBColor(15)
            .Text23(24).BackColor = QBColor(15)
            .Text23(23).BackColor = QBColor(15)
            .Text23(22).BackColor = QBColor(15)

            If TOFreqFlag = 1 Then
                .Text23(25).Enabled = True
                .Text23(24).Enabled = True
                .Text23(23).Enabled = True
                .Text23(22).Enabled = True
                .Label21(28).Enabled = True
                .Label21(27).Enabled = True
                .Label21(26).Enabled = True
                .Label21(25).Enabled = True
            Else
                .Text23(25).Enabled = False
                .Text23(24).Enabled = False
                .Text23(23).Enabled = False
                .Text23(22).Enabled = False
                .Label21(28).Enabled = False
                .Label21(27).Enabled = False
                .Label21(26).Enabled = False
                .Label21(25).Enabled = False
            End If

        End If

        '  TOTreeType As Integer, TOPerms As Integer, TOWinLen As Integer, TOStepSize As Integer, TOSmooth As Integer
        'Public TOHigh As Double, MatAverage As Double, TOPower As Double, TOPValCOff As Double, TOTvTs As Double, TOFreqFlag As Double, TOFreqA As Double, TOFreqC As Double, TOFreqG As Double, TOFreqT As Double
        .Text1(23) = TOSmooth
        .Text1(26) = TORndNum
        .Text1(25) = TOPower
        .Text1(27) = TOPerms
        .Text1(24) = TOPValCOff

        If TOTreeType = 0 Then
            .Label1(34) = "Construct only LS trees"
        Else
            .Label1(34) = "Construct NJ and LS trees"
        End If

        xTOWinLen = TOWinLen
        xToTreeType = TOTreeType
        xToFreqFlag = TOFreqFlag
        xTOTsTv = TOTvTs
        xTOModel = TOModel
        'Tree Options
        .Text1(33) = TBSReps
        .Text1(34) = TRndSeed
        .Text1(32) = TTVRat


        
        

        If TModel = 0 Then
            .Label1(44).Caption = "Jukes and Cantor, 1969"
            .Frame21(6).Enabled = False
            .Text1(32).Text = "0.5"
            .Text1(32).Enabled = False
            .Text1(32).BackColor = Form1.BackColor
            .Label1(43).Enabled = False
            .Text1(31).Enabled = False
            .Text1(31).BackColor = Form1.BackColor
            .Text1(31).ForeColor = QBColor(8)
            .Label1(42).Enabled = False
            .Label21(35).Enabled = False
            .Command28(26).Enabled = False
            .Text23(30).BackColor = Form1.BackColor
            .Text23(31).BackColor = Form1.BackColor
            .Text23(32).BackColor = Form1.BackColor
            .Text23(33).BackColor = Form1.BackColor
            .Text23(30).Enabled = False
            .Text23(31).Enabled = False
            .Text23(32).Enabled = False
            .Text23(33).Enabled = False
            .Label21(36).Enabled = False
            .Label21(37).Enabled = False
            .Label21(38).Enabled = False
            .Label21(39).Enabled = False
        ElseIf TModel = 1 Then
            .Label1(44).Caption = "Kimura, 1980"
            .Text1(32).Text = TTVRat
            .Text1(32).Enabled = True
            .Text1(32).BackColor = QBColor(15)
            .Label1(43).Enabled = True
            .Text1(31).Enabled = False
            .Text1(31).BackColor = Form1.BackColor
            .Text1(31).ForeColor = QBColor(8)
            .Label1(42).Enabled = False
            .Frame21(6).Enabled = False
            .Label21(35).Enabled = False
            .Command28(26).Enabled = False
            .Text23(30).BackColor = Form1.BackColor
            .Text23(31).BackColor = Form1.BackColor
            .Text23(32).BackColor = Form1.BackColor
            .Text23(33).BackColor = Form1.BackColor
            .Text23(30).Enabled = False
            .Text23(31).Enabled = False
            .Text23(32).Enabled = False
            .Text23(33).Enabled = False
            .Label21(36).Enabled = False
            .Label21(37).Enabled = False
            .Label21(38).Enabled = False
            .Label21(39).Enabled = False
        ElseIf TModel = 2 Then
            .Label1(44).Caption = "Jin  and  Nei, 1990"
            .Text1(32).Text = TTVRat
            .Text1(32).Enabled = True
            .Text1(32).BackColor = QBColor(15)
            .Label1(43).Enabled = True
            .Text1(31).Enabled = True
            .Text1(31).BackColor = QBColor(15)
            .Text1(31).ForeColor = QBColor(0)
            .Label1(42).Enabled = True
            .Frame21(6).Enabled = False
            .Label21(35).Enabled = False
            .Command28(26).Enabled = False
            .Text23(30).BackColor = Form1.BackColor
            .Text23(31).BackColor = Form1.BackColor
            .Text23(32).BackColor = Form1.BackColor
            .Text23(33).BackColor = Form1.BackColor
            .Text23(30).Enabled = False
            .Text23(31).Enabled = False
            .Text23(32).Enabled = False
            .Text23(33).Enabled = False
            .Label21(36).Enabled = False
            .Label21(37).Enabled = False
            .Label21(38).Enabled = False
            .Label21(39).Enabled = False
        ElseIf TModel = 3 Then
            .Label1(44).Caption = "Felsenstein, 1984"
            .Text1(32).Text = TTVRat
            .Text1(32).Enabled = True
            .Text1(32).BackColor = QBColor(15)
            .Label1(43).Enabled = True
            .Text1(31).Enabled = False
            .Text1(31).BackColor = Form1.BackColor
            .Text1(31).ForeColor = QBColor(8)
            .Label1(42).Enabled = False
            .Frame21(6).Enabled = True
            .Label21(35).Enabled = True
            .Command28(26).Enabled = True
            .Text23(30).BackColor = QBColor(15)
            .Text23(31).BackColor = QBColor(15)
            .Text23(32).BackColor = QBColor(15)
            .Text23(33).BackColor = QBColor(15)

            If TBaseFreqFlag = 1 Then
                .Text23(30).BackColor = QBColor(15)
                .Text23(31).BackColor = QBColor(15)
                .Text23(32).BackColor = QBColor(15)
                .Text23(33).BackColor = QBColor(15)
                .Text23(30).Enabled = True
                .Text23(31).Enabled = True
                .Text23(32).Enabled = True
                .Text23(33).Enabled = True
                .Label21(36).Enabled = True
                .Label21(37).Enabled = True
                .Label21(38).Enabled = True
                .Label21(39).Enabled = True
            Else
                .Text23(30).BackColor = Form1.BackColor
                .Text23(31).BackColor = Form1.BackColor
                .Text23(32).BackColor = Form1.BackColor
                .Text23(33).BackColor = Form1.BackColor
                .Text23(30).Enabled = False
                .Text23(31).Enabled = False
                .Text23(32).Enabled = False
                .Text23(33).Enabled = False
                .Label21(36).Enabled = False
                .Label21(37).Enabled = False
                .Label21(38).Enabled = False
                .Label21(39).Enabled = False
            End If

        End If

        .Text1(31) = TCoeffVar
        '     TBaseFreqFlag
        .Text23(31) = TAfreq
        .Text23(32) = TCFreq
        .Text23(33) = TGFreq
        .Text23(30) = TTFreq
        .Text1(35) = TPower
        ' Write #1, TNegBLFlag,

        If TBaseFreqFlag = 0 Then
            .Label21(35).Caption = "Estimate from alignment"
        Else
            .Label21(35).Caption = "User defined"
        End If

        If TNegBLFlag = 0 Then
            .Label1(49) = "Negative branch lengths not allowed"
        Else
            .Label1(49) = "Negative branch lengths allowed"
        End If

        If TSubRepsFlag = 0 Then
            .Label1(50) = "Do not do subreplicates"
        Else
            .Label1(50) = "Do subreplicates"
        End If

        If TGRFlag = 0 Then
            .Label1(51) = "Do not do global rearrangements"
        Else
            .Label1(51) = "Do global rearrangements"
        End If

        If TRndIOrderFlag = 0 Then
            .Label1(52) = "Do not randomise input order"
        Else
            .Label1(52) = "Randomise input order"
        End If
        
        
        xTModel = TModel
        xTBaseFreqFlag = TBaseFreqFlag
        
        xTPModel = TPModel
        xTPBPFEstimate = TPBPFEstimate
        .Text1(17) = TPTVRat
        .Text1(18) = TPInvSites
        .Text23(34) = TPGamma
        .Text23(42) = TPAlpha
         xModelTestFlag = ModelTestFlag
        
        Call SetMLModel
       
        If ModelTestFlag = 0 Then
            .Label1(69) = "User specified model"
            
        Else
            .Label1(69) = "Automatic model selection"
            .Label1(24).Caption = "Jukes and Cantor, 1969"
            .Frame21(6).Enabled = False
            
            'Disable Transition transversion ratio option
            .Text1(17).Text = "0.5"
            .Text1(17).Enabled = False
            .Text1(17).BackColor = Form1.BackColor
            .Label1(59).Enabled = False
            
            'Disable base frequency estimate option
            .Label21(46).Enabled = False
            .Command28(44).Enabled = False
            
            
            
            .Text23(30).BackColor = Form1.BackColor
            .Text23(31).BackColor = Form1.BackColor
            .Text23(32).BackColor = Form1.BackColor
            .Text23(33).BackColor = Form1.BackColor
            .Text23(30).Enabled = False
            .Text23(31).Enabled = False
            .Text23(32).Enabled = False
            .Text23(33).Enabled = False
            .Label21(36).Enabled = False
            .Label21(37).Enabled = False
            .Label21(38).Enabled = False
            .Label21(39).Enabled = False
        End If
        
        
        
        xTBModel = TBModel
        xTBGamma = TBGamma
        Call SetTBModel
        Call SetTBGamma
        
        
        .Text1(20) = TBGammaCats
        .Text1(39) = TBNGens
        .Text1(40) = TBNChains
        .Text1(38) = TBSampFreq
        .Text1(41) = TBTemp
        .Text1(42) = TBSwapFreq
        .Text1(43) = TBSwapNum
        
        
        
        
        
        xTNegBLFlag = TNegBLFlag
        xTSubRepsFlag = TSubRepsFlag
        xTGRFlag = TGRFlag
        xTRndIOrderFlag = TRndIOrderFlag
        
        If OptionsFlag = 0 Then
            Dim OnIndex As Byte
            If .TabStrip1.SelectedItem.Index = 1 Then
                OnIndex = 1
            ElseIf .TabStrip1.SelectedItem.Index = 2 Then
                OnIndex = 2
            ElseIf .TabStrip1.SelectedItem.Index = 3 Then
                OnIndex = 3
            ElseIf .TabStrip1.SelectedItem.Index = 4 Then
                OnIndex = 4
            ElseIf .TabStrip1.SelectedItem.Index = 5 Then
                OnIndex = 5
            ElseIf .TabStrip1.SelectedItem.Index = 6 Then
                OnIndex = 6
            ElseIf .TabStrip1.SelectedItem.Index = 7 Then
                OnIndex = 7
            ElseIf .TabStrip1.SelectedItem.Index = 8 Then
                OnIndex = 13
            ElseIf .TabStrip1.SelectedItem.Index = 9 Then
                OnIndex = 14
            ElseIf .TabStrip1.SelectedItem.Index = 10 Then
                OnIndex = 8
            ElseIf .TabStrip1.SelectedItem.Index = 11 Then
                OnIndex = 11 '9
            ElseIf .TabStrip1.SelectedItem.Index = 12 Then
                OnIndex = 10
            ElseIf .TabStrip1.SelectedItem.Index = 13 Then
                OnIndex = 9 '11
            ElseIf .TabStrip1.SelectedItem.Index = 14 Then
                OnIndex = 12
            ElseIf .TabStrip1.SelectedItem.Index = 15 Then
                OnIndex = 1
            ElseIf .TabStrip1.SelectedItem.Index = 16 Then
                OnIndex = 1
            End If
            For X = 0 To 13

                If X = OnIndex - 1 Then
                    .Frame2(X).Visible = True
                    'Else
                    '    Frame2(X).Visible = False
                Else
                    .Frame2(X).Visible = False
                End If

            Next 'X

        End If

    End With
End Sub
Public Sub VarRecRates(Nextno, StrainSeq() As String, GCFlag, GCTractLen, BlockPen, StartRho, MCMCUpdates, FreqCo, FreqCoMD)

'MCMCUpdates = 100000
Dim FN As String, SubSample As Long, DoSeqNo As Long, TSeq() As String, HS As Long, DS() As Byte, NewNo As Long, Discard As Long

ReDim Preserve RefCol(0)

Call UnModNextno
Form1.Picture2.Enabled = False

ManFlag = 20

'Exit Sub
    Form1.SSPanel1.Caption = "Loading Likelihood Lookup"
    
    Rnd (-BSRndNumSeed)
    Discard = MCMCUpdates / 10
    ODir = CurDir
    ChDir App.Path
    XX = CurDir
    'Check for the available likelihood lookup files and choose one
        MaxSize = 99 'set maxsize according to the max size of the lookup file
        
    
    
    
    'Get sequences to analyse
    'currently can handle only 100 sequences
    
    If Nextno > MaxSize Then
        'if over 100 sequences in the alignment then a random sample of 100 is looked at
        DoSeqNo = MaxSize
        ReDim DS(Nextno)
        Randomize
        ReDim TSeq(MaxSize)
        For X = 0 To MaxSize
            Do
                NewNo = Int((Nextno * Rnd) + 1)
                
                If DS(NewNo) = 0 Then
                    DS(NewNo) = 1
                    Exit Do
                End If
            Loop
            TSeq(X) = StrainSeq(NewNo)
        Next X
    Else
        DoSeqNo = Nextno
        ReDim TSeq(Nextno)
        
        For X = 0 To Nextno
            TSeq(X) = StrainSeq(X)
        Next X
        
    End If
    
    
    
    If CircularFlag = 1 Then
        HS = CInt(Len(StrainSeq(0)) / 2)
        For X = 0 To DoSeqNo
            TSeq(X) = right$(TSeq(X), HS) + TSeq(X) + left$(TSeq(X), HS)
            XX = Len(TSeq(X))
        Next X
    End If
    
    
    'Make sequence file
    Close #1
    Open "seqs.fas" For Output As #1
    Print #1, Str(DoSeqNo + 1) + " " + Str(Len(TSeq(0))) + " 1"
    For X = 0 To DoSeqNo
        Print #1, ">S" + Trim(Str(X))
        Print #1, TSeq(X)
    Next X
    Close #1
    'vrflag = 1
    Dim FC As Long
    Form1.ProgressBar1.Value = 2
    
    'Make batch file to run "infiles"
    On Error Resume Next
    Kill "of"
    On Error GoTo 0
    Open "convert.bat" For Output As #1
    Print #1, "convert < infile > of"
    Close #1
    'Make infile for convert
    'first do it without the frequency cutoff
    Open "infile" For Output As #1
    Print #1, "seqs.fas"
    Print #1, "0"
    Print #1, "2"
    Print #1, "0"
    Print #1, Str(FreqCoMD)
    Print #1, "1"
    Print #1, "0"
    Close #1
    'XX = CurDir
    'run convert
    Call ShellAndClose("convert.bat", 0)
    
    Open "of" For Binary Access Read As #1
    Dim InString As String, TS As String
    
    InString = String(LOF(1), " ")
    Get #1, , InString
    Close #1
    Pos = InStr(1, InString, "Summary of output data", vbBinaryCompare)
    For X = 0 To 5
        LastPos = Pos + 1
        Pos = InStr(LastPos, InString, "=", vbBinaryCompare)
        TS = Mid$(InString, Pos + 1, 10)
        TS = Trim(TS)
        VarRho(4 + X) = Val(TS)
        If Pos = 0 Then Exit For
    Next X
    
    
    
    
    
    Form1.ProgressBar1.Value = 5
    Dim MSF As String
    'given theta choose the most suitable likelihood file
    MSF = "LF"
    
    If VRFlag = 0 Then
        Open "lkgen.bat" For Output As #1
        Print #1, "lkgen < infile"
        Close #1
        FN = "LF" + Trim(Str(DoSeqNo + 1))
        Open FN For Append As #1
        
        If LOF(1) = 0 Then
            Close #1
            z = DoSeqNo + 2
            Do
                ChDir App.Path
                FN = "LF" + Trim(Str(z))
                Open FN For Append As #1
                FC = LOF(1)
                Close #1
                If FC > 0 Then
                    Exit Do
                ElseIf z < 100 Then
                    Kill FN
                End If
                z = z + 1
                XX = CurDir
                If z > 100 Then
                    Response = MsgBox("No likelihood lookup file was found.  Go to the LDHat web site, download a likelihood lookup file, rename it 'LF100' and copy it to " + App.Path, 0, "RDP Warning")
                
                    Exit Sub
                    
                End If
            Loop
            'make likelihood file
            'make infile for lkgen
            Call LoadLF(FN, DoSeqNo)
            If ErrorFlag = 1 Then
                ErrorFlag = 0
                Close #1
                Exit Sub
            End If
        Else
            Close #1
        End If
        
    End If
    
    'load it
    
    
    
    'now do it with the frequency cutoff
    On Error Resume Next
    Kill "of"
    On Error GoTo 0
    Open "convert.bat" For Output As #1
    Print #1, "convert < infile > of"
    Close #1
    'Make infile for convert
    Open "infile" For Output As #1
    Print #1, "seqs.fas"
    Print #1, "0"
    Print #1, "2"
    Print #1, Str(FreqCo)
    Print #1, Str(FreqCoMD)
    Print #1, "1"
    Print #1, "0"
    Close #1
    'XX = CurDir
    'run convert
    Call ShellAndClose("convert.bat", 0)
    
    
    
    If GCFlag = 1 Then
        Open "locs" For Binary As #1
        Dim WholeFile As String
        WholeFile = String(LOF(1), " ")
        Get #1, , WholeFile
        Close #1
        
        Pos = InStr(1, WholeFile, "L", vbBinaryCompare)
        If Pos > 0 Then
            Mid$(WholeFile, Pos, 1) = "C"
        End If
        Open "locs" For Output As 1
        Print #1, WholeFile
        Close #1
    End If
    
    If VRFlag = 0 Then
        'Make batch file to run "infiles"
        Open "interval.bat" For Output As #1
        Print #1, "interval < infile"
        Close #1
        
        'Make "infile for interval"
        Open "infile" For Output As #1
        
        Print #1, "sites"
        Print #1, "locs"
        If GCFlag = 1 Then
            Print #1, Trim(Str(GCTractLen))
        End If
        Print #1, FN
        Print #1, "0"
        Print #1, Str(StartRho)
        Print #1, Str(BlockPen)
        Print #1, Str(MCMCUpdates)
        
        If MCMCUpdates / 1000 > 2000 Then
            SubSample = 2000
        Else
            SubSample = CInt(MCMCUpdates / 1000)
            If SubSample < 1 Then SubSample = 1
        End If
        
        Discard = CLng(Discard / SubSample)
        Print #1, Str(SubSample)
        
        Close #1
        
        Form1.SSPanel1.Caption = "0 of " + Trim(Str(MCMCUpdates)) + " MCMC Updates Completed"
        'Now execute and wait till its finished
        Dim NSites As Long
        
        Open "locs" For Input As #1
        Input #1, NSites
        Close #1
        
        TargetFileSize = NSites
        'set size of rates to 0
        Open "rates.txt" For Output As #1
        Close #1
        
        
        TargetFileSize = (TargetFileSize * 9) * (MCMCUpdates / SubSample)
        
        BatIndex = 8
        FullSize = MCMCUpdates
        On Error Resume Next
        Kill "rates.txt"
        On Error GoTo 0
        Dim Prob As Byte
        Prob = 0
        Do
        
            Call ShellAndClose("interval.bat", 0)
            Open "rates.txt" For Append As #1
            If LOF(1) > 10 Then
                
                Close #1
                Exit Do
            Else
                If Prob = 0 Then
                'there is probably a problem with the likelihood file
                    Prob = 1
                    On Error Resume Next
                    Kill FN
                    On Error GoTo 0
                    Call LoadLF(FN, DoSeqNo)
                    If ErrorFlag = 1 Then
                        ErrorFlag = 0
                        Form1.SSPanel1.Caption = ""
                        Form1.ProgressBar1 = 0
                        
                        Exit Sub
                    End If
                Else
                    MsgBox "There has been a problem executing the Interval component of LDHat"
                    Form1.SSPanel1.Caption = ""
                    Form1.ProgressBar1.Value = 0
                    Screen.MousePointer = 0
                    Close #1
                    Exit Sub
                End If
            End If
            Close #1
        Loop
        'Make batch file to run "infiles"
        Open "stat.bat" For Output As #1
        Print #1, "stat < infile"
        Close #1
        Open "stat.bat" For Output As #1
        Print #1, "stat rates.txt " + Trim(Str(Discard))
        Close #1
        
        Call ShellAndClose("stat.bat", 0)
        X = X
    Else
        Open "locs" For Input As #1
        Input #1, NSites
        Close #1
    End If
        'Close #1
        Open "Res.txt" For Input As #1
        'read in results
        Dim Crap As String, CP As Double
        Line Input #1, Crap
        Line Input #1, Crap
        
        Dim PltVals() As Double
        ReDim PltVals(2, NSites + 1)
        For X = 1 To NSites
            Input #1, CP
            Input #1, PltVals(0, X - 1)
            Input #1, CP
            Input #1, CP
            Input #1, PltVals(1, X - 1)
            Input #1, PltVals(2, X - 1)
            
        Next X
        Close #1

    'pltvals 0,1 is where the actual data begins
    
    Dim DPX() As Long, PDX() As Long, Max As Double
    ReDim DPX(Len(TSeq(0))), PDX(Len(TSeq(0)))
    Open "locs" For Input As #1
    Line Input #1, Crap
    For X = 1 To NSites
        Input #1, DPX(X)
        If DPX(X) <= Len(TSeq(0)) Then
            PDX(DPX(X)) = X
        End If
    Next X
    Close #1
    Max = 0
    If CircularFlag = 1 Then
        SP = HS
        EP = Len(TSeq(0)) - HS
    Else
        SP = 0
        EP = Len(TSeq(0)) - 1
    End If
    Dim GraphPlt() As Double, PltPos() As Double, LPP As Double, TotSite As Long
    ReDim GraphPlt(2, Len(TSeq(0))), PltPos(Len(TSeq(0))), XDiffpos(Len(TSeq(0)))
    Y = 0
    
    TotSite = 0
    
    
    zz = 0
    For X = 1 To NSites
        
        If DPX(X) >= SP And X > 1 Then
            
            
            If PltVals(2, X) > Max Then Max = PltVals(2, X)
            GraphPlt(0, Y) = PltVals(0, X - 1)
            GraphPlt(1, Y) = PltVals(1, X - 1)
            GraphPlt(2, Y) = PltVals(2, X - 1)
            
            PltPos(Y) = (DPX(X) - SP)
            'If Y > 0 Then
                LPP = (DPX(X - 1) - SP)
            'End If
            TotSite = TotSite + (PltPos(Y) - LPP)
            VarRho(0) = VarRho(0) + (PltPos(Y) - LPP) * PltVals(0, X - 1)
            VarRho(1) = VarRho(1) + (PltPos(Y) - LPP) * PltVals(1, X - 1)
            VarRho(2) = VarRho(2) + (PltPos(Y) - LPP) * PltVals(2, X - 1)
            If (DPX(X - 1) - SP) > 0 Then
                zz = zz + 1
                XDiffpos(zz) = (DPX(X - 1) - SP)
            
            End If
            
            If X = NSites And X = 12345678 Then
                PltPos(Y) = PltPos(Y - 1) + (PltPos(Y) - LPP) / 2
                If PltPos(Y) > Len(StrainSeq(0)) Then
                    PltPos(Y) = PltPos(Y) - Len(StrainSeq(0))
                    'Move everything over 1
                    Dim TmpPos() As Double
                    ReDim TmpPos(3)
                    TmpPos(0) = PltPos(Y)
                    TmpPos(1) = GraphPlt(0, Y - 1)
                    TmpPos(2) = GraphPlt(1, Y - 1)
                    TmpPos(3) = GraphPlt(2, Y - 1)
                    For z = NSites - 1 To 1 Step -1
                        PltPos(z) = PltPos(z - 1)
                        GraphPlt(0, z) = GraphPlt(0, z - 1)
                        GraphPlt(1, z) = GraphPlt(1, z - 1)
                        GraphPlt(2, z) = GraphPlt(2, z - 1)
                    Next z
                    PltPos(0) = TmpPos(0)
                    GraphPlt(0, 0) = TmpPos(1)
                    GraphPlt(0, 1) = TmpPos(2)
                    GraphPlt(0, 2) = TmpPos(3)
                End If
            Else
                
                'If Y > 0 Then
                    PltPos(Y) = LPP + (PltPos(Y) - LPP) / 2
                'Else
                   ' PltPos(Y) = PltPos(Y) / 2
                'End If
                X = X
            End If
            If DPX(X) > EP Then
                If PltPos(Y) > Len(StrainSeq(0)) Then
                    PltPos(Y) = Len(StrainSeq(0))
                End If
                Y = Y + 1
                Exit For
            End If
            Y = Y + 1
        Else
            X = X
        End If
    Next X
    NSites = Y
For X = 0 To 2
    VarRho(X) = VarRho(X) / TotSite
Next X
VarRho(3) = TotSite
If CircularFlag = 1 Then
    VarRho(4) = VarRho(4) / 2
    VarRho(6) = VarRho(6) / 2
End If
Dim Pict As Long, PntAPI As POINTAPI, YScaleFactor As Double
YScaleFactor = 0.85
PicHeight = Form1.Picture7.Height * YScaleFactor
XFactor = ((Form1.Picture7.Width - 40) / Len(StrainSeq(0)))
LenXoverSeq = NSites
Call DoAxes(0, Len(StrainSeq(0)), -1, Max, 0, 1, "Rho(4Ner) per bp")

XFactor = ((Form1.Picture7.Width - 40) / Len(StrainSeq(0)))
Pict = Form1.Picture7.hDC

'get everything into standard save/copy format
ReDim GPrint(1, NSites * 2 + 2), GPrintCol(1), GPrintPos(1, NSites * 2 + 2)

ReDim GVarPos(0, NSites)
For X = 1 To NSites
    GVarPos(0, X) = XDiffpos(X)
Next X
ReDim GCritval(10)
GLegend = "Rho(4Ner) per bp"
GPrintLen = NSites * 2 + 2 'how many points to plot
GPrintCol(0) = 0 'line is black
GPrintCol(1) = RGB(128, 128, 128) 'line is grey
GPrintNum = 1 'two lines
GPrintType = 0 'a normal line plot
GPrintMin(0) = 0 'bottom val
GPrintMin(1) = Max 'upper val

For X = 0 To NSites - 1
    
    GPrint(0, X) = GraphPlt(0, X)
    GPrint(0, NSites * 2 - X) = GraphPlt(0, X)
    GPrintPos(0, X) = PltPos(X)
    GPrintPos(0, NSites * 2 - X) = PltPos(X)
    'GVarPos(0,X) = PltPos(X)
    
Next X

GBlockNum = -1



For X = 0 To NSites - 1
    
    GPrint(1, X) = GraphPlt(1, X)
    GPrint(1, NSites * 2 - X) = GraphPlt(2, X)
    GPrintPos(1, X) = PltPos(X)
    GPrintPos(1, NSites * 2 - X) = PltPos(X)
    
Next X

For X = 0 To NSites - 1
    For Y = 0 To 1
    
        If GPrintPos(Y, X) < 1 Then GPrintPos(Y, X) = 1
    Next Y
    
Next X
GExtraTNum = 1
ReDim GExtraText(GExtraTNum)
GExtraText(0) = "Recombination rate"
GExtraText(1) = "95% Confidence interval"

GPrintPos(1, GPrintLen - 1) = GPrintPos(1, 0)
GPrintPos(0, GPrintLen - 1) = GPrintPos(0, 0)
GPrintPos(1, GPrintLen) = GPrintPos(1, 0)
GPrintPos(0, GPrintLen) = GPrintPos(0, 0)
GPrint(0, GPrintLen - 1) = GPrint(0, 0)
GPrint(1, GPrintLen - 1) = GPrint(1, GPrintLen - 2)
GPrint(1, GPrintLen) = GPrint(1, 0)
Form1.Picture7.AutoRedraw = True
Form1.Frame17.Visible = False
If X = X Then
    Form1.Picture7.DrawStyle = 0
    Form1.Picture7.ForeColor = RGB(150, 150, 150)
    'MoveToEx Pict, 30 + PltPos(0) * XFactor, PicHeight - (15 + GraphPlt(1, 0) / Max * (PicHeight - 35)), PntAPI
    'LineTo Pict, 30 + PltPos(X) * XFactor, (PicHeight - (15 + (GraphPlt(2, 0) / Max) * (PicHeight - 35)))
    For X = 0 To NSites - 2
        
        If PltPos(X) > 0 And PltPos(X) <= Len(StrainSeq(0)) Then
            ST = CLng(((GraphPlt(2, X)) / Max) * (PicHeight - 35) - ((GraphPlt(1, X)) / Max) * (PicHeight - 35))
            EN = CLng(((GraphPlt(2, X + 1)) / Max) * (PicHeight - 35) - ((GraphPlt(1, X + 1)) / Max) * (PicHeight - 35))
            
            If ST = 0 Then ST = 1
            If EN = 0 Then EN = 1
            If ST > EN Then
                For Y = 0 To ST
                    MoveToEx Pict, 30 + PltPos(X) * XFactor, PicHeight - ((15 + GraphPlt(1, X) / Max * (PicHeight - 35))) - ST * (Y / ST), PntAPI
                    LineTo Pict, 30 + PltPos(X + 1) * XFactor, (PicHeight - ((15 + (GraphPlt(1, X + 1) / Max) * (PicHeight - 35)))) - EN * (Y / ST)
                Next Y
            Else
                
                For Y = 0 To EN
                
                    MoveToEx Pict, 30 + PltPos(X) * XFactor, PicHeight - ((15 + GraphPlt(1, X) / Max * (PicHeight - 35))) - ST * (Y / EN), PntAPI
                    LineTo Pict, 30 + PltPos(X + 1) * XFactor, (PicHeight - ((15 + (GraphPlt(1, X + 1) / Max) * (PicHeight - 35)))) - EN * (Y / EN)
                    X = X
                Next Y
            End If
        End If
        
    Next 'X
    
    
    
    Form1.Picture7.ForeColor = 0
    Y = 0
    For z = 0 To NSites - 1
        If PltPos(z) > 0 Then
            MoveToEx Pict, 30 + PltPos(z) * XFactor, PicHeight - (15 + GraphPlt(Y, 0) / Max * (PicHeight - 35)), PntAPI
            Exit For
        End If
    Next z
    For X = 1 To NSites - 1
        
        LineTo Pict, 30 + PltPos(X) * XFactor, (PicHeight - (15 + (GraphPlt(Y, X) / Max) * (PicHeight - 35)))
    Next 'X
    
Else
    For Y = 0 To 2
        If Y > 0 Then
            Form1.Picture7.DrawStyle = 2
        Else
            Form1.Picture7.DrawStyle = 0
        End If
        MoveToEx Pict, 30 + PltPos(0) * XFactor, PicHeight - (15 + GraphPlt(Y, 0) / Max * (PicHeight - 35)), PntAPI
        
        For X = 1 To NSites - 1
        
            LineTo Pict, 30 + PltPos(X) * XFactor, (PicHeight - (15 + (GraphPlt(Y, X) / Max) * (PicHeight - 35)))
        Next 'X
    Next Y
End If
ManFlag = 20
'Call DoLegend
 VRFlag = 1
Form1.Picture7.DrawStyle = 0
Form1.ProgressBar1.Value = 0
Form1.SSPanel1.Caption = ""
ChDir ODir
Form1.Picture2.Enabled = True
Form1.Picture21.PaintPicture Form1.Picture7.Image, Form1.Picture7.left, Form1.Picture7.Top + 5
End Sub
Public Sub NumToString(Num, SD, OutString As String, Exp As String, ExtraX As String)
Dim TString As String, TempNum As Double, SMod As Byte
TString = Trim(Str(Num))
If Num < 0 Then SMod = 3 Else SMod = 2
Exp = ""
If Num < 1 And Num > -1 Then 'ie num is a decimal
    'check if it is a very small number
    Pos = InStr(1, TString, "E", vbBinaryCompare)
    If Pos > 0 Then 'ie its a very small number
        Exp = Mid$(TString, Pos + 1, Len(TString) - Pos)
        If Pos > SD + SMod Then
            TempNum = Val(left(TString, Pos - 1))
            TempNum = TempNum * (10 ^ (SD))
            TempNum = CLng(TempNum)
            TempNum = TempNum / (10 ^ (SD))
            TString = Trim(Str(TempNum))
            Pos = InStr(1, TString, "E", vbBinaryCompare)
            If Pos = 0 Then
                TString = left$(TString, SD + SMod)
            Else
                TString = left(TString, Pos - 1)
                
            End If
        Else
            TString = left(TString, Pos - 1)
        End If
    Else
        'check for a large number of decimal places
        TempNum = Num * 10 ^ SD
        If TempNum < 1 And TempNum > -1 Then
            TempNum = Num
            TempNum = CLng(-Log10(TempNum))
            tempnum2 = TempNum
            TempNum = Num * 10 ^ TempNum
            Exp = Trim(Str(-tempnum2))
            TempNum = TempNum * (10 ^ (SD))
            TempNum = CLng(TempNum)
            TempNum = TempNum / (10 ^ (SD))
            TString = Trim(Str(TempNum))
            X = X
        Else
            Exp = ""
            TString = "0" + TString
        End If
    End If
Else
    TempNum = Num
    TempNum = TempNum * (10 ^ SD)
    TempNum = CLng(TempNum)
    TempNum = TempNum / (10 ^ SD)
    TString = Trim(Str(TempNum))
End If
If Exp = "" Then ExtraX = "" Else ExtraX = " X 10"
OutString = TString

End Sub
Public Sub MakeAPos(APos() As Long, BPos() As Long)
Dim TypeSeq As Long
ReDim APos(Len(StrainSeq(0))), BPos(Len(StrainSeq(0)))
    TypeSeq = TypeSeqNumber
    If TypeSeq > PermNextNo Or TypeSeq < 0 Then TypeSeq = 0
    
    For X = 1 To Len(StrainSeq(0))
        'If X = Len(StrainSeq(0)) - 10 Then
        '    X = X
        'End If
        APos(X) = X - SeqSpaces(X, TypeSeq)
        BPos(X - SeqSpaces(X, TypeSeq)) = X
    Next X
End Sub
Public Sub RMinWin()

Dim ST As Long, EN As Long, Win As Long, STX As Long, ENX As Long, RMP() As Long, MaxBP As Double, UScore As Long, A As Long, B As Long
Win = 200
ReDim RMP(Len(StrainSeq(0)))
Y = 0
z = 0
If DoneMatX(6) = 0 Then
    Call RMinMat(GCFlag, GCTractLen, FreqCo, FreqCoMD, 6)
    CurMatrixFlag = 6
    Call DoMatCap
    If (MatZoom(CurMatrixFlag) + 0.5) * Form1.Picture26.ScaleHeight - Form1.Picture26.ScaleHeight <= 32000 Then
        Form1.Command40.Enabled = True
    End If
    
    Form1.Command39.Enabled = True
End If
For X = Win / 2 + 1 To Len(StrainSeq(0)) - Win / 2
    If X = 2735 Then
        X = X
    End If
    
    ST = X - Win / 2
    EN = X + Win / 2 - 1
    
    Do While VarsitesLD(Y) <= ST
        Y = Y + 1
        If Y > UBound(VarsitesLD) Then Exit For
    Loop
    STX = Y 'VarsitesLD(Y)
    Do While VarsitesLD(z) <= EN
        z = z + 1
        If z > Mat567Len - 1 Then z = Mat567Len - 1: Exit Do
        If z > UBound(VarsitesLD) Then Exit For
    Loop
    z = z - 1
    ENX = z 'VarsitesLD(Z)
    If STX > UBound(MatrixRMin, 1) Then STX = UBound(MatrixRMin, 1)
    If ENX > UBound(MatrixRMin, 2) Then ENX = UBound(MatrixRMin, 2)
 
    RMP(X) = MatrixRMin(CLng((STX / Mat567Len) * UBound(MatrixRMin, 1)), CLng((ENX / Mat567Len) * UBound(MatrixRMin, 1)))
    X = X
    'RMP(X) = MatrixRMin(STX, ENX)
 'If X = X Then
 '       If ST = 1465 Then
 '           X = X
 '       End If
 '
 '       If MatrixRMin(STX, ENX) > 0 Then
 '           X = X
 '       End If
 '   Else
 '       UScore = 0
 '       For A = STX To ENX
 '           For B = STX To ENX
 '               If UScore < MatrixRMin(A, B) Then
 '                   UScore = MatrixRMin(A, B)
 '               End If
 '           Next B
 '
 '       Next A
 '       RMP(X) = UScore
 '   End If
    
Next X

MaxBP = 0
For X = 0 To Len(StrainSeq(0))
    If MaxBP < RMP(X) Then
        MaxBP = RMP(X)
    End If

Next X

If MaxBP = 0 Then Exit Sub

Form1.Picture7.AutoRedraw = True
Form1.Picture7.Picture = LoadPicture()
YScaleFactor = 0.85
PicHeight = Form1.Picture7.Height * YScaleFactor

Dim PntAPI As POINTAPI
Dim Pict As Long

Dim BPos() As Long, TypeSeq
Pict = Form1.Picture7.hDC




Call MakeAPos(APos(), BPos())
XFactor = ((Form1.Picture7.Width - 40) / APos(Len(StrainSeq(0))))
Call DoAxes(1, Len(StrainSeq(0)), -1, MaxBP, 0, 0, "Breakpoints per" + Str(Win) + " nt window")

Form1.Picture7.DrawWidth = 3


Form1.Picture7.ForeColor = RGB(128, 128, 128)
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + APos((Win / 2 + 1)) * XFactor, PicHeight - 15 - RMP(Win / 2 + 1) / MaxBP * (PicHeight - 35), PntAPI

For X = Win / 2 + 1 To Len(StrainSeq(0)) - Win / 2 - 1
    LineTo Pict, 30 + APos(X) * XFactor, PicHeight - 15 - RMP(X) / MaxBP * (PicHeight - 35)
Next 'X



If RelX > 0 Or RelY > 0 Then

    Call Highlight

End If
Form1.Picture7.DrawWidth = 1
Form1.Picture7.ForeColor = 0
MoveToEx Pict, 30 + APos((Win / 2 + 1)) * XFactor, PicHeight - 15 - RMP(Win / 2 + 1) / MaxBP * (PicHeight - 35), PntAPI

For X = Win / 2 + 1 To Len(StrainSeq(0)) - Win / 2 - 1
    LineTo Pict, 30 + APos(X) * XFactor, PicHeight - 15 - RMP(X) / MaxBP * (PicHeight - 35)
Next 'X


Form1.Picture7.Refresh
End Sub
Public Sub RMinMat(GCFlag, GCTractLen, FreqCo, FreqCoMD, MatFlag)

Dim Pict As Long, PntAPI As POINTAPI, YScaleFactor As Double, MaxN As Double
Dim FN As String, SubSample As Long, DoSeqNo As Long, TSeq() As String, HS As Long, DS() As Byte, NewNo As Long, Discard As Long
If DoneMatX(5) = 0 Then
    ReDim Preserve RefCol(0)
    Form1.Picture2.Enabled = False
    ManFlag = 20
    Form1.SSPanel1.Caption = "Loading Likelihood Lookup"
    
    Rnd (-BSRndNumSeed)
    Discard = MCMCUpdates / 10
    ODir = CurDir
    ChDir App.Path
    'Get sequences to analyse
    'currently can handle only 100 sequences
    
    If Nextno > 99 Then
        'if over 100 sequences in the alignment then a random sample of 100 is looked at
        DoSeqNo = 99
        ReDim DS(Nextno)
        For X = 0 To 99
            Do
                NewNo = Int((Nextno * Rnd) + 1)
                
                If DS(NewNo) = 0 Then
                    DS(NewNo) = 1
                    Exit Do
                End If
            Loop
            TSeq(X) = StrainSeq(X)
        Next X
    Else
        DoSeqNo = Nextno
        ReDim TSeq(Nextno)
        
        For X = 0 To Nextno
            TSeq(X) = StrainSeq(X)
        Next X
        
    End If
    
    'Make sequence file
    Open "seqs.fas" For Output As #1
    Print #1, Str(DoSeqNo + 1) + " " + Str(Len(TSeq(0))) + " 1"
    For X = 0 To DoSeqNo
        Print #1, ">S" + Trim(Str(X))
        Print #1, TSeq(X)
    Next X
    Close #1
    
    Open "lkgen.bat" For Output As #1
    Print #1, "lkgen < infile"
    Close #1
    
    FN = "LF" + Trim(Str(DoSeqNo + 1))
    Open FN For Append As #1
    If LOF(1) = 0 Then
        Close #1
        'make likelihood file
        'make infile for lkgen
        Open "infile" For Output As #1
        Print #1, "lf100"
        Print #1, Trim(Str(DoSeqNo + 1))
        Close #1
        
        Call ShellAndClose("lkgen.bat", 0)
        On Error Resume Next
        Kill FN
        On Error GoTo 0
        Name "new_lk.txt" As FN
    Else
        Close #1
    End If
    Form1.ProgressBar1.Value = 5
    'Make batch file to run "infiles"
    Open "convert.bat" For Output As #1
    Print #1, "convert < infile" ' > outfilec"
    Close #1
    'Make infile for convert
    Open "infile" For Output As #1
    Print #1, "seqs.fas"
    Print #1, "0"
    Print #1, "2"
    Print #1, Str(FreqCo)
    Print #1, Str(FreqCoMD)
    Print #1, "1"
    Print #1, "0"
    Close #1
    
    'run convert
    Call ShellAndClose("convert.bat", 0)
    
    If GCFlag = 1 Then
        Open "locs" For Binary As #1
        Dim WholeFile As String
        WholeFile = String(LOF(1), " ")
        Get #1, , WholeFile
        Close #1
    
        Pos = InStr(1, WholeFile, "L", vbBinaryCompare)
        Mid$(WholeFile, Pos, 1) = "C"
        Open "locs" For Output As 1
        Print #1, WholeFile
        Close #1
    
    End If
   ' Exit Sub
    'Make batch file to run "infiles"
    Open "pairwise.bat" For Output As #1
    Print #1, "pairwise < infile > outfilep"
    Close #1
    
    'Make "infile for pairwise"
    Open "infile" For Output As #1
    
    Print #1, "sites"
    Print #1, "locs"
    If GCFlag = 1 Then
        Print #1, Trim(Str(GCTractLen))
    End If
    Print #1, "1"
    Print #1, FN
    Print #1, "0"
    
    Print #1, "0"
    
    Print #1, "0"
    
    Print #1, "2"
    Print #1, "0"
    Print #1, "0"
    Print #1, "0"
    Close #1
    
    'Form1.SSPanel1.Caption = "0 of " + Trim(Str(MCMCUpdates)) + " MCMC Updates Completed"
    'Now execute and wait till its finished
    
    
    Call ShellAndClose("pairwise.bat", 0)
    
    
    
    Open "locs" For Input As #1
    Input #1, NSites
    Close #1
    Dim tMat() As Double
    If NSites - 1 < 1000 Then
        RSize = 1000
    Else
        RSize = NSites - 1
    End If
    ReDim MatrixLD(RSize + 1, RSize + 1)
    ReDim MatrixRMin(RSize + 1, RSize + 1)
    ReDim MatrixRMinD(RSize + 1, RSize + 1)
    ReDim tMat(RSize + 1, RSize + 1)
    
    'get site positions
    Dim DPX() As Long, PDX() As Long, Max As Double
    ReDim DPX(Len(TSeq(0))), PDX(Len(TSeq(0)))
    Open "locs" For Input As #1
    Line Input #1, Crap
    ReDim VarsitesLD(Len(StrainSeq(0)))
    For X = 1 To NSites
        Input #1, DPX(X)
        PDX(DPX(X)) = X
        VarsitesLD(X) = DPX(X)
    Next X
    
    Close #1
    
    Max = 0
    SP = 1
    EP = Len(TSeq(0))
    
    Dim GraphPlt() As Double, PltPos() As Double
    ReDim GraphPlt(2, Len(TSeq(0))), PltPos(Len(TSeq(0))), XDiffpos(Len(TSeq(0)))
    Y = 0
    For X = 1 To NSites
        If DPX(X) >= SP And DPX(X) <= EP Then
            
            
            PltPos(Y) = DPX(X) - (SP - 1)
            XDiffpos(Y) = PltPos(Y)
            Y = Y + 1
        Else
            X = X
        End If
    Next X
    'open LDmatrix
    Open "fit.txt" For Input As #1
        Line Input #1, Crap
        Line Input #1, Crap
        Line Input #1, Crap
        Line Input #1, Crap
        Y = NSites - 2
        For X = 0 To NSites - 1
            Input #1, CP
            For z = 0 To Y
                Input #1, tMat((z + ((NSites - 1) - Y) + 1), X)
                tMat(X, z + ((NSites - 2) - Y) + 1) = MatrixLD(z + ((NSites - 2) - Y) + 1, X)
                
            Next z
            Y = Y - 1
        Next X
    Close #1
    Dim ScaleFact As Double
    ScaleFact = RSize / NSites
    If ScaleFact > 1 Then
        ScaleFact = ScaleFact
    End If
    For Y = 0 To RSize
        B = CLng(Y / ScaleFact)
        For X = Y + 1 To RSize
            A = CLng(X / ScaleFact)
            
            MatrixLD(Y, X) = tMat(A, B)
            MatrixLD(X, Y) = tMat(A, B)
            
        Next X
    Next Y
    ReDim tMat(RSize + 1, RSize + 1)
    Open "rmin.txt" For Input As #1
        Line Input #1, Crap
        Line Input #1, Crap
        Line Input #1, Crap
        Line Input #1, Crap
        Line Input #1, Crap
        For X = 0 To NSites - 2
            Input #1, CP
            
            For z = 0 To NSites - 1
                If z <> X Then
                    Input #1, tMat(z, X)
                    X = X
                End If
            Next z
            
        Next X
    Close #1
    
    For Y = 0 To RSize
        B = CLng(Y / ScaleFact)
        For X = Y + 1 To RSize
            A = CLng(X / ScaleFact)
            
            MatrixRMin(Y, X) = tMat(A, B)
            MatrixRMin(X, Y) = tMat(A, B)
            MatrixRMinD(Y, X) = tMat(B, A)
            MatrixRMinD(X, Y) = tMat(B, A)
            X = X
        Next X
    Next Y
    
    MaxN = FindMaxN(RSize, MatrixLD(0, 0))
    MatBound(5) = MaxN
    
    MaxN = 0
    For X = 0 To RSize
        For Y = X + 1 To RSize
            If MatrixLD(X, Y) < MaxN Then MaxN = MatrixLD(X, Y)
        Next Y
    Next X
    If Abs(MaxN) > MatBound(5) Then MatBound(5) = Abs(MaxN)
    
    MaxN = FindMaxN(RSize, MatrixRMin(0, 0))
    MatBound(6) = MaxN
    'For X = 0 To RSize
    '    For Y = 0 To RSize
    '        If MatrixRMin(X, Y) > 0 Then
    '            X = X
    '        End If
    '    Next Y
    'Next X
    For X = 0 To RSize
        For Y = 0 To RSize
            If MatrixRMinD(X, Y) > 0.5 Then
                MatrixRMinD(X, Y) = 0.5 'Then
             End If
        Next Y
    Next X
    MaxN = FindMaxN(RSize, MatrixRMinD(0, 0))
    
    'MaxN = 0
    'For X = 0 To NSites
    '    For Y = 0 To NSites
    '        If MatrixRMinD(X, Y) > MaxN Then
    '            MaxN = MatrixRMinD(X, Y)
    '        End If
    '    Next Y
    'Next X
    
    MatBound(7) = MaxN
    DoneMatX(5) = 1
    Mat567Len = NSites
Else

    NSites = Mat567Len
    
    Max = 0
    
    
    If NSites - 1 < 1000 Then
        RSize = 1000
    Else
        RSize = NSites - 1
    End If
   
    
End If

'RSize = NSites - 1
'MatBound(MatFlag) = 0.1
MaxN = MatBound(MatFlag)

Form1.Picture26.Picture = LoadPicture()
    
    
Dim PosS(1) As Long, PosE(1) As Long, DistD As Long
XAddj = (Form1.Picture26.ScaleHeight) / RSize
DistD = RSize / MatZoom(MatFlag)
PosS(0) = MatCoord(MatFlag, 0)
PosE(0) = PosS(0) + DistD
PosS(1) = MatCoord(MatFlag, 1)
PosE(1) = PosS(1) + DistD
Dim Limit As Long
If MatFlag = 5 Then
    Limit = UBound(MatrixLD, 1)
ElseIf MatFlag = 6 Then
    Limit = UBound(MatrixRMin, 1)
ElseIf MatFlag = 7 Then
    Limit = UBound(MatrixRMinD, 1)
End If
    
If PosE(1) > (Limit - 1) - 1 Then PosE(1) = (Limit - 1) - 1
If PosE(0) > (Limit - 1) - 1 Then PosE(0) = (Limit - 1) - 1
       
If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
        
If MatFlag = 6 Then
    If MaxN > 2 Then
        MaxN = CLng((MaxN / 4) + 0.5) * 4
    Else
        MaxN = 2
    End If
End If

Form1.Picture26.Picture = LoadPicture()
Form1.Picture26.ScaleMode = 3
DistD = RSize / MatZoom(MatFlag)
XAddj = (Form1.Picture26.ScaleHeight) / DistD
     
If MatFlag = 5 Then
    Call DrawmatsVB(-MaxN, PosE(), PosS(), SX, SY, XAddj, MatrixLD(), HeatMap(), CurScale, MaxN)
    If DontDoKey = 0 Then
        Call DoKey(1, MaxN, -MaxN, MatFlag, "Marginal likelihood ratio", CurScale)
    End If
ElseIf MatFlag = 6 Then
    Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixRMin(), HeatMap(), CurScale, MaxN)
    If DontDoKey = 0 Then
        Call DoKey(1, MaxN, 0, MatFlag, "RMin", CurScale)
    End If
    
ElseIf MatFlag = 7 Then
    Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixRMinD(), HeatMap(), CurScale, MaxN)
    If DontDoKey = 0 Then
        Call DoKey(1, MaxN, 0, MatFlag, "RMin/D", CurScale)
    End If
End If
    Form1.Picture26.Refresh
    Form1.Picture26.Enabled = True
    Form1.SSPanel6(2).Enabled = True
    Form1.Picture29.Enabled = True
    Form1.ProgressBar1 = 0
    Form1.SSPanel1 = ""
End Sub
Public Sub CentreBP(BE, EN, XPosDiff() As Long, XDiffpos() As Long, BWarn, EWarn, LSS, LenXoverSeq)


If XPosDiff(BE) - 1 > 0 Then
    BE = BE - CLng(((BE - XDiffpos(XPosDiff(BE) - 1)) / 2) - 0.1)
Else
    BE = BE - CLng(((BE + Len(StrainSeq(0)) - XDiffpos(LenXoverSeq)) / 2) - 0.1)
End If

If BE = 0 Then
    BE = 1
ElseIf BE < 1 Then
    If CircularFlag = 0 Then 'seqs linear
        BE = 1
    Else
        BE = Len(StrainSeq(0)) + BE
    End If
End If

z = BE
If SEventNumber > 0 Then
    If MissingData(z, Seq1) = 1 Or MissingData(z, Seq2) = 1 Or MissingData(z, Seq3) = 1 Then
        BWarn = 1
        Do
            z = z + 1
            
            If z > Len(StrainSeq(0)) Then
                If CircularFlag = 0 Then
                    z = 1
                    Exit Do
                Else
                    z = 1
                End If
            End If
            If MissingData(z, Seq1) = 0 And MissingData(z, Seq2) = 0 And MissingData(z, Seq3) = 0 Then
               
                Exit Do
            End If
        Loop
    End If
End If
BE = z
If CircularFlag = 0 Then
    If XPosDiff(BE) < LSS Then
        BWarn = 1
    End If
End If
XPosDiff(Len(StrainSeq(0))) = LenXoverSeq
'XX = XDiffpos(LenXoverSeq)
'EN = 850

If XPosDiff(EN) + 1 <= LenXoverSeq Then
    EN = EN + CLng(((XDiffpos(XPosDiff(EN) + 1) - EN) / 2) - 0.1)
Else
    EN = EN + CLng(((XDiffpos(1) + (Len(StrainSeq(0)) - EN)) / 2) - 0.1)
End If
If EN > Len(StrainSeq(0)) Then
                                     
    If CircularFlag = 0 Then 'seqs linear
        EN = Len(StrainSeq(0))
        
    Else
        EN = EN - Len(StrainSeq(0))
    End If
End If
z = EN
If SEventNumber > 0 Then
    If MissingData(z, Seq1) = 1 Or MissingData(z, Seq2) = 1 Or MissingData(z, Seq3) = 1 Then
        EWarn = 1
        Do
            z = z - 1
            
            If z < 1 Then
                If CircularFlag = 0 Then
                    z = Len(StrainSeq(0))
                    Exit Do
                Else
                    z = Len(StrainSeq(0))
                End If
            End If
            If MissingData(z, Seq1) = 0 And MissingData(z, Seq2) = 0 And MissingData(z, Seq3) = 0 Then
               
                Exit Do
            End If
        Loop
    End If
    EN = z
End If
If CircularFlag = 0 Then
    If XPosDiff(EN) > LenXoverSeq - LSS Then
        EWarn = 1
    End If
End If
End Sub
Public Sub DrawCompatMat()
Dim PSize As Long, Pict As Long
Dim X As Long, Y As Long, RelaventSites() As Long, tCnt As Byte, VarNum As Long, BinMat() As Byte
Dim CompSeq As Long, tCMat() As Byte

If DoneMatX(0) = 0 Then
    DoneMatX(0) = 1
    ReDim RelaventSites(100, Len(StrainSeq(0)))
    ReDim VarsitesCM(Len(StrainSeq(0)))
    For Y = 1 To Len(StrainSeq(0))
        For X = 0 To Nextno
            RelaventSites(SeqNum(Y, X), Y) = RelaventSites(SeqNum(Y, X), Y) + 1
        Next X
    Next Y
    VarNum = 0
    For Y = 1 To Len(StrainSeq(0))
        tCnt = 0
        If RelaventSites(66, Y) > 1 Then
            tCnt = tCnt + 1
        End If
        If RelaventSites(68, Y) > 1 Then
            tCnt = tCnt + 1
        End If
        If RelaventSites(72, Y) > 1 Then
            tCnt = tCnt + 1
        End If
        If RelaventSites(85, Y) > 1 Then
            tCnt = tCnt + 1
        End If
        If tCnt = 2 Then
            VarsitesCM(VarNum) = Y
            VarNum = VarNum + 1
        End If
    Next Y
    VarNum = VarNum - 1
    ReDim BinMat(VarNum, Nextno)
    For X = 0 To VarNum
        For Y = 0 To Nextno
            If SeqNum(VarsitesCM(X), Y) <> 46 Then
                CompSeq = Y
                Exit For
            End If
        Next Y
        For Y = 0 To Nextno
            If SeqNum(VarsitesCM(X), Y) <> 46 Then
                If SeqNum(VarsitesCM(X), Y) <> SeqNum(VarsitesCM(X), CompSeq) Then
                    BinMat(X, Y) = 1
                Else
                    BinMat(X, Y) = 0
                End If
            Else
                BinMat(X, Y) = 2
            End If
        Next Y
    Next X
    ReDim tCMat(1, 1)
    ReDim MatrixC(VarNum, VarNum)
    Dim NTS As Long, NDone As Long
    NTS = (VarNum + 1) * (VarNum) / 2
    For X = 0 To VarNum
        SS = GetTickCount
        For Y = X + 1 To VarNum
            tCMat(0, 0) = 0: tCMat(0, 1) = 0: tCMat(1, 0) = 0: tCMat(1, 1) = 0
            NDone = NDone + 1
            GoOn = 0
            For z = 0 To Nextno
                A = BinMat(X, z)
                B = BinMat(Y, z)
                'Exit Sub
                If A < 2 And B < 2 Then
                    If tCMat(A, B) = 0 Then
                        tCMat(A, B) = 1
                        tCnt = tCMat(0, 0) + tCMat(0, 1) + tCMat(1, 0) + tCMat(1, 1)
                        If tCnt = 4 Then GoOn = 1: Exit For
                    End If
                End If
            Next z
            'tCnt = tCMat(0, 0) + tCMat(0, 1) + tCMat(1, 0) + tCMat(1, 1)
            If GoOn = 1 Then
                MatrixC(X, Y) = 1
                MatrixC(Y, X) = 1
            End If
            
        Next Y
        EE = GetTickCount
            TT = EE - SS
            X = X
        nt = Abs(GetTickCount)
        If Abs(nt - LT) > 500 Then
            LT = nt
            Form1.ProgressBar1 = 10 + (NDone / NTS) * 70
            Form1.SSPanel1.Caption = Str(NDone) + " of" + Str(NTS) + " site pairs screened"
            Form1.Refresh
        End If
    Next X
    MatBound(0) = VarNum
Else
    VarNum = MatBound(0)
End If
    
SS = GetTickCount

    
    Dim cAddj As Long, MatPic() As Double, PosS(1) As Long, PosE(1) As Long, DistD As Long
    
    
    DistD = VarNum / MatZoom(0)
    
    Form1.Picture26.ScaleMode = 3
    XAddj = (Form1.Picture26.ScaleHeight) / DistD
    cAddj = (Int((1 / XAddj) + 1)) ^ 2
    cAddj = Int(255 / cAddj) '(255 * (XAddj ^ 2)) - 1
    If cAddj > 255 Then cAddj = 255
    If XAddj <= 1 Then
        PSize = XAddj * VarNum
    Else
        PSize = VarNum
    End If
    If XAddj > 1 Then
        Span = Int(XAddj + 1)
    Else
        Span = 1
    End If
    Dim SpT As Long, SpE As Long
    SpT = CLng(-(Span / 2)) + 1
    SpE = CLng((Span / 2) - 0.00001)
    If SpT > SpE Then SpT = SpE
    
    'if spt=spe
    PosS(0) = MatCoord(0, 0)
    PosE(0) = PosS(0) + DistD
    PosS(1) = MatCoord(0, 1)
    PosE(1) = PosS(1) + DistD
    If PosE(1) > (UBound(MatrixC, 1) - 1) - 1 Then PosE(1) = (UBound(MatrixC, 1) - 1) - 1
    If PosE(0) > (UBound(MatrixC, 1) - 1) - 1 Then PosE(0) = (UBound(MatrixC, 1) - 1) - 1
       
    If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
    If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
    ReDim MatPic(PSize + 1, PSize + 1)
    'If XAddj < 1 Then
    SS = GetTickCount
    If XAddj <= 1 Then
        If SpT <> SpE Then
            For Y = SY To PosE(1)
                For X = SX To PosE(0)
                
                    If MatrixC(X, Y) = 1 Then
                        For A = SpT To SpE
                            If CLng(X * XAddj) + A >= 0 And CLng(X * XAddj) + A <= PSize Then
                                For B = SpT To SpE
                                    If CLng(Y * XAddj) + B >= 0 And CLng(Y * XAddj) + B <= PSize Then
                                        If MatPic(CLng(X * XAddj + 0.001) + A, CLng(Y * XAddj + 0.001) + B) + cAddj <= 255 Then
                                            MatPic(CLng(X * XAddj + 0.001) + A, CLng(Y * XAddj + 0.001) + B) = MatPic(CLng(X * XAddj + 0.001) + A, CLng(Y * XAddj + 0.001) + B) + cAddj
                                            'MatPic(CLng(Y * XAddj) + B, CLng(X * XAddj) + A) = MatPic(CLng(X * XAddj) + A, CLng(Y * XAddj) + B)
                                        End If
                                    End If
                                Next B
                            End If
                        Next A
                    End If
                Next X
            Next Y
        
        
        
        Else
            
            For Y = SY To PosE(1)
                YP = CLng(Y * XAddj + 0.001)
                For X = SX To PosE(0)
                
                    If MatrixC(X, Y) = 1 Then
                        XP = CLng(X * XAddj + 0.001)
                        'For A = SpT To SpE
                            'If CLng(X * XAddj) + A >= 0 And CLng(X * XAddj) + A <= PSize Then
                                'For B = SpT To SpE
                                    'If CLng(Y * XAddj) + B >= 0 And CLng(Y * XAddj) + B <= PSize Then
                                        'If MatPic(XP, YP) + cAddj <= 255 Then
                                            MatPic(XP, YP) = MatPic(XP, YP) + cAddj
                                            'MatPic(CLng(Y * XAddj) + B, CLng(X * XAddj) + A) = MatPic(CLng(X * XAddj) + A, CLng(Y * XAddj) + B)
                                        'End If
                                    'End If
                                'Next B
                            'End If
                        'Next A
                    End If
                Next X
            Next Y
        End If
    Else
        For X = 0 To VarNum
            For Y = 0 To VarNum
                If MatrixC(X, Y) = 1 Then
                    MatPic(X, Y) = 255
                End If
            Next Y
        Next X
    End If
    EE = GetTickCount
    TT = EE - SS
'Else
    SS = GetTickCount
    If X = 1234 Then 'this is not working for matpic
        MaxN = FindMaxN(PSize, MatPic(0, 0))
        
    Else
        MaxN = 0
        For X = 0 To PSize
            For Y = 0 To PSize
                If MatPic(X, Y) > MaxN Then
                    MaxN = MatPic(X, Y)
                End If
            Next Y
        Next X
        
    End If
    MaxN = MaxN / cAddj
    If MaxN = 0 Then MaxN = 1
    EE = GetTickCount
    TT = EE - SS
'Else
'    For X = 0 To VarNum
'        For Y = 0 To VarNum
'            If MatrixC(X, Y) = 1 Then
'                MatPic(CLng(X * XAddj), CLng(Y * XAddj)) = MatPic(CLng(X * XAddj), CLng(Y * XAddj)) + cAddj
'            End If
'        Next Y
'    Next X
'End If
Form1.Picture26.Picture = LoadPicture()
Form1.Picture26.ScaleMode = 3
Form1.Picture26.AutoRedraw = True


'If XAddj > 1 Then
'    Form1.Picture26.DrawWidth = CInt(XAddj + 1)
'End If
Pict = Form1.Picture26.hDC
Dim PntAPI As POINTAPI
'If XAddj < 1 Then xaddjx = 1 Else xaddjx = XAddj
SS = GetTickCount
Dim PosEx(1) As Long, PosSx(1) As Long

For X = 0 To 1
    PosEx(X) = CLng(PosE(X) * XAddj + 0.001)
    PosSx(X) = CLng(PosS(X) * XAddj + 0.001)
Next X
If X = X Then
    If XAddj <= 1 Then
        Call DrawmatsVB(0, PosEx(), PosSx(), SX * XAddj, SY * XAddj, 1, MatPic(), HeatMap, 1, MaxN * cAddj)
    Else
        Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatPic(), HeatMap, 1, 255)
    End If
Else
    
    For X = (SX * XAddj) To PosEx(0)
        For Y = SY * XAddj To PosEx(1) 'Step 20\
        
            'If MatPic(X, Y) > 0 Or X = X Then
                z = X + 1
                RM = MatPic(X, Y)
                Do While z <= PosEx(0)
                    If MatPic(z, Y) <> RM Then Exit Do
                    z = z + 1
                Loop
                z = z - 1
                If z = X Then
                    SetPixelV Pict, (X - PosSx(0)), Y - PosSx(1), HeatMap(1, (MatPic(X, Y)) / 255 * 1020)
                Else
                    YP = Y - PosSx(1)
                    Form1.Picture26.ForeColor = HeatMap(1, (MatPic(X, Y)) / 255 * 1020)
                    Dummy = MoveToEx(Pict, (X - PosSx(0)), YP, PntAPI)
                    Dummy = LineTo(Pict, (z - PosSx(0)) + 1, YP)
                    X = X
                    'Form1.Picture26.Refresh
                End If
                X = X
                'SetPixelV Pict, Y - PosS(0) - PosS(1), X, HeatMap(1, (MatPic(Y, X)) / 255 * 1020)
            'End If
        Next Y
        'Form1.Picture26.Refresh
    Next X
End If
 EE = GetTickCount
 TT = EE - SS
    X = X
X = X
If DontDoKey = 0 Or X = X Then
    Call DoKey(0, MaxN, 0, 0, "Number of incompatible sites", 1)
End If
Form1.Picture26.Refresh
Form1.SSPanel1.Caption = ""
Form1.ProgressBar1 = 0
End Sub
Public Sub DoColourScale()
Dim PHeight As Long
Form3.Picture1.ScaleMode = 3
Form3.Picture1.AutoRedraw = True
PHeight = Form3.Picture1.Height
XAddj = Form3.Picture1.ScaleWidth / 1020

For X = 0 To 1020
    Form3.Picture1.Line (X * XAddj, 0)-(X * XAddj, PHeight), HeatMap(CurScale, X)
Next X
Form3.Picture1.Refresh

End Sub
Public Sub DoAlnAddj(SNNextNo, AdjustAlign() As Long)

Dim TDist() As Double, WinSeq(1) As Long, MaxDist As Double, SeqsIn() As Integer, SeqsOut() As Integer, MaxLen As Long
ReDim TDist(Nextno)
'do a profile aignment of sequences in two separate RDP files

'First work out which sequences have the smallest distance to all others
For X = 0 To SNNextNo - 1
    For Y = 0 To SNNextNo - 1
        TDist(X) = TDist(X) + Distance(X, Y)
    Next Y
Next X
For X = SNNextNo To Nextno
    For Y = SNNextNo To Nextno
        TDist(X) = TDist(X) + Distance(X, Y)
    Next Y
Next X
MaxDist = 0
For X = 0 To SNNextNo - 1
    If MaxDist < TDist(X) Then
        MaxDist = TDist(X)
        WinSeq(0) = X
    End If
Next X

MaxDist = 0
For X = SNNextNo To Nextno
    If MaxDist < TDist(X) Then
        MaxDist = TDist(X)
        WinSeq(1) = X
    End If
Next X

MaxLen = Len(StrainSeq(0)) * 2
ReDim SeqsIn(MaxLen, 1), SeqsOut(MaxLen, 1)

For z = 0 To 1
    Y = 0
    For X = 1 To Len(StrainSeq(0))
        'If SeqNum(X, WinSeq(Z)) <> 46 Then
            Y = Y + 1
            SeqsIn(Y, z) = SeqNum(X, (WinSeq(z)))
        'End If
    Next X
Next z
Form1.SSPanel1.Caption = "Doing profile alignment"
Call doAlignmentSh(SeqsIn(), SeqsOut(), 1, MaxLen, 10, 5, 1)
Dim tSeqNum() As Integer, NonG(1, 1) As Long, NonGSeq()
NonG(0, 0) = 0: NonG(1, 0) = 0: NonG(0, 1) = 0: NonG(1, 1) = 0


ReDim tSeqNum(MaxLen, Nextno), NonGSeq(MaxLen, 1, 1)
For z = 0 To 1
    For X = 1 To MaxLen
        If X <= UBound(SeqNum, 1) Then
            If SeqNum(X, WinSeq(z)) <> 46 Then
                NonG(z, 0) = NonG(z, 0) + 1
                NonGSeq(NonG(z, 0), z, 0) = X
            End If
        End If
        If SeqsOut(X, z) <> 46 And SeqsOut(X, z) <> 0 Then
            NonG(z, 1) = NonG(z, 1) + 1
            NonGSeq(NonG(z, 1), z, 1) = X
        End If
    Next X
Next z
Dim AddGaps() As Long, SP(1), EP(1)
ReDim AddGaps(MaxLen, 1)

SP(0) = 0: SP(1) = SNNextNo
EP(0) = SNNextNo - 1: EP(1) = Nextno

ReDim AdjustAlign(Len(StrainSeq(0)) * 2, 1)


For z = 0 To 1
    For X = 1 To NonG(z, 1)
        If NonGSeq(X + 1, z, 0) > 0 Then
            A = NonGSeq(X + 1, z, 0) - NonGSeq(X, z, 0) - 1
        Else
            A = 0
        End If
        
        For B = 0 To A
            AdjustAlign(NonGSeq(X, z, 0) + B, z) = NonGSeq(X, z, 1) + B
            'tSeqNum(NonGSeq(X, Z, 1) + B, Y) = SeqNum(NonGSeq(X, Z, 0) + B, Y)
        Next B
        For Y = SP(z) To EP(z)
            
            For B = 0 To A
                tSeqNum(NonGSeq(X, z, 1) + B, Y) = SeqNum(NonGSeq(X, z, 0) + B, Y)
            Next B
        Next Y
    Next X
    
Next z




Dim EndSeq As Long

For z = MaxLen To 1 Step -1
    If tSeqNum(z, 0) <> 0 Then
        For X = z To 1 Step -1
            GoOn = 1
            For Y = 0 To Nextno
                If tSeqNum(z, Y) <> 46 Then
                    GoOn = 0
                    Exit For
                End If
            Next Y
            If GoOn = 0 Then
                EndSeq = X
                Exit For
            End If
        Next X
        Exit For
    End If
Next z

'For Y = 1 To MaxLen
'    If SeqsOut(Y, 0) = 0 Then
'        EndSeq = Y - 1
'        Exit For
'    End If
'Next Y


ReDim SeqNum(EndSeq, Nextno)
For X = 0 To Nextno
    For Y = 1 To EndSeq
        If tSeqNum(Y, X) <> 0 Then
            SeqNum(Y, X) = tSeqNum(Y, X)
        Else
            SeqNum(Y, X) = 46
        End If
    Next Y
Next X

ReDim StrainSeq(Nextno)
Dim TSeq As String
For X = 0 To Nextno
    TSeq = ""
    For Y = 1 To EndSeq
        TSeq = TSeq + Chr(SeqNum(Y, X) - 1)
    Next Y
    StrainSeq(X) = TSeq
Next X
XX = UBound(AdjustAlign, 1)
For z = 0 To 1
    For X = SP(z) To EP(z)
        For Y = 1 To CurrentXOver(X)
            XOverList(X, Y).Beginning = AdjustAlign(XOverList(X, Y).Beginning, z)
            XOverList(X, Y).Ending = AdjustAlign(XOverList(X, Y).Ending, z)
        Next Y
        For Y = 1 To BCurrentXoverMa(X)
            BestXOListMa(X, Y).Beginning = AdjustAlign(BestXOListMa(X, Y).Beginning, z)
            BestXOListMa(X, Y).Ending = AdjustAlign(BestXOListMa(X, Y).Ending, z)
        Next Y
        For Y = 1 To BCurrentXoverMi(X)
            BestXOListMi(X, Y).Beginning = AdjustAlign(BestXOListMi(X, Y).Beginning, z)
            BestXOListMi(X, Y).Ending = AdjustAlign(BestXOListMi(X, Y).Ending, z)
        Next Y
    Next X
Next z

End Sub
Public Sub RecombMapPermsB(DN As Long, APos() As Long, BPos() As Long, Excl() As Byte, CV() As Double, Win As Long, PermNum As Long, RNDSeed As Long)
    Dim SPS As Long, EPS As Long, A As Long, X As Long, MaxPos As Long, MinPos As Long, OffsetX As Long, PPos As Long, Target(1) As Long, MaxV(1) As Double, MaxP(1) As Long, MaxVals() As Long, NewStart As Long, s As Long, NS As Long, NE As Long, LSSeq As Long, RecSize As Long, MaxS As Long, Map() As Integer, Size As Long, d As Long, P1 As Long, P2 As Long, RecMapSmooth() As Double
    Dim tSeqNum() As Integer, MaxCycleNo As Long, CycleNo As Long
    Dim BestP() As Double, BreakPos() As Long, BreakNum As Long, OKProg(10) As Byte
    Dim xxWin As Double
    ReDim RecMap(Len(StrainSeq(0)))
    ReDim RecMapSmooth(Len(StrainSeq(0)))
    If DoScans(0, 0) = 1 Then OKProg(0) = 1
    If DoScans(0, 1) = 1 Then OKProg(1) = 1
    If DoScans(0, 2) = 1 Then OKProg(2) = 1
    If DoScans(0, 3) = 1 Then OKProg(3) = 1
    If DoScans(0, 4) = 1 Then OKProg(4) = 1
    If DoScans(0, 5) = 1 Then OKProg(5) = 1
    'Dim MinS As Long
    
    If CircularFlag = 0 Then
        MaxPos = Len(StrainSeq(0)) - Win / 2
        MinPos = Win / 2
    Else
        MaxPos = Len(StrainSeq(0))
        MinPos = 1
    End If
    
    If LongWindedFlag = 0 Then
        ENumb = Eventnumber
    Else
        ENumb = SEventNumber
    End If
    B = 0
    C = 0
    ReDim Map(Len(StrainSeq(0)), PermNum), MaxVals(1, PermNum)
    Rnd (-RNDSeed)
    SS = GetTickCount
    ReDim XDiffpos(Len(StrainSeq(0)) + 200), XPosDiff(Len(StrainSeq(0)) + 200)
    For A = 1 To PermNum
        ReDim tSeqNum(Len(StrainSeq(0)), Nextno)
        If CircularFlag = 0 Then
            For X = 0 To Nextno
            '    For Y = 1 To Len(StrainSeq(0))
                   tSeqNum(1, X) = 1
                   tSeqNum(Len(StrainSeq(0)), X) = 1
            '    Next Y
            Next X
        End If
        zzz = 0
        ggg = 0
        
        
        For X = 1 To ENumb
            
            
            If Excl(X) = 1 Then
                zzz = zzz + 1
                If (BestEvent(X, 0) > 0 Or BestEvent(X, 1) > 0) Then
                    
                    CNum = 0
                    sen = SuperEventList(XOverList(BestEvent(X, 0), BestEvent(X, 1)).Eventnumber)
                    
                    
                    
                    If X = X Then
                        d = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Daughter
                        P1 = XOverList(BestEvent(X, 0), BestEvent(X, 1)).MajorP
                        P2 = XOverList(BestEvent(X, 0), BestEvent(X, 1)).MinorP
                        ST = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Beginning
                        EN = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Ending
                        
                        'Make variable sites array
                        
                        
                        LSSeq = MakeSubSeqPerm(Len(StrainSeq(0)), d, P1, P2, PermSeqNum(0, 0), XDiffpos(0), XPosDiff(0))
                        
                        If EN > ST Then
                            RecSize = XPosDiff(EN) - XPosDiff(ST) + 1
                        Else
                            RecSize = LSSeq - XPosDiff(ST) + XPosDiff(EN)
                        End If
                        If CircularFlag = 0 Then
                            MaxS = LSSeq - RecSize
                        Else
                             MaxS = LSSeq
                        End If
                        
                        'Make new ending and start for this event
                    
                        'For A = 1 To PermNum
                            MaxCycleNo = LSSeq * 5
                            CycleNo = 0
                            Do
                                NewStart = Int((MaxS * Rnd) + 1)
                               
                                GoOn = CheckBPOL(Len(StrainSeq(0)), d, LSSeq, NewStart, RecSize, tSeqNum(0, 0), NS, NE, XDiffpos(0))
                                If GoOn = 1 Then Exit Do
                                CycleNo = CycleNo + 1
                                If CycleNo > MaxCycleNo Then Exit Do
                                
                            Loop
                            OffsetX = NS - ST
                            If XOverList(BestEvent(X, 0), BestEvent(X, 1)).SBPFlag <> 3 Then
                                
                                    For d = 0 To 1
                                        GoOn = 0
                                        If d = 0 Then
                                            
                                            If XOverList(BestEvent(X, 0), BestEvent(X, 1)).SBPFlag <> 1 Then
                                                s = APos(NS): GoOn = 1
                                               
                                            End If
                                        Else
                                            If XOverList(BestEvent(X, 0), BestEvent(X, 1)).SBPFlag <> 2 Then
                                                s = APos(NE): GoOn = 1
                                            End If
                                        End If
                                    
                                        If GoOn = 1 Then
                                            
                                            Dummy = AddToMap(A, s, Win, Len(StrainSeq(0)), APos(0), Map(0, 0))
                                            
                                        End If
                                    Next d
                                
                            End If
                            
                         'Next A
                    
                   
                        
        
                    End If
                    If X = X Then
                        'ReDim BestP(1, NextNo), BreakPos(1, 1, NextNo)
                        
                        'mark BPs for all the other events
                        'For Y = 1 To SEventNumber
                            SPS = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Beginning + OffsetX
                            EPS = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Ending + OffsetX
                            If SPS < 1 Then SPS = Len(StrainSeq(0)) + SPS
                            If SPS > Len(StrainSeq(0)) Then SPS = SPS - Len(StrainSeq(0))
                            If EPS < 1 Then EPS = Len(StrainSeq(0)) + EPS
                            If EPS > Len(StrainSeq(0)) Then EPS = EPS - Len(StrainSeq(0))
                            For z = 0 To PermNextNo
                                If Daught(X, z) > 0 Then
                                    tSeqNum(SPS, z) = 1
                                    tSeqNum(EPS, z) = 1
                                End If
                            Next z
                        'Next Y
                    End If
                    
                    
                End If
            End If
        Next X
        ET = GetTickCount
        If Abs(ET - LT) > 500 Then
            LT = ET
            Form1.SSPanel1.Caption = Str(A) + " of" + Str(PermNum) + " permutations completed"
            Form1.ProgressBar1.Value = (A / PermNum) * 80
            If A / 10 = CLng(A / 10) Then
                Form1.Refresh
            Else
                Form1.SSPanel1.Refresh
            End If
        End If
        xxx = zzz
    Next A
    'get maxes
   EE = GetTickCount
   TT = EE - SS
   '16.141
   X = X
   'Exit Sub
    'DN = 1
    SS = GetTickCount
    If X = X Then
        Dummy = FindMaxMapVal(DN, Len(StrainSeq(0)), PermNum, MaxVals(0, 0), Map(0, 0))
    Else
        For X = 1 To PermNum
            For Y = DN To Len(StrainSeq(0)) - DN + 1
                If MaxVals(0, X) < Map(Y, X) Then MaxVals(0, X) = Map(Y, X)
               
            Next Y
            X = X '9,8,9,9,8,8,8,12,7,8,8,10,10
        Next X
    End If
    EE = GetTickCount
    TT = EE - SS
    '3.194
    'XX = MaxVals(0, 5) '16,16,16,19,20
    SS = GetTickCount
    
    'sort through the top 5% and get the 99% and 95% highest scores
    
    Target(1) = CLng(0.05 * PermNum) + 1
    Target(0) = CLng(0.01 * PermNum) + 1
    For z = 0 To 1
        If Target(z) < 1 Then Target(z) = 1
    Next z
    For z = 1 To Target(1)
        MaxV(0) = 0
        MaxV(1) = 0
        For X = 0 To PermNum
            If MaxV(0) < MaxVals(0, X) Then
                MaxV(0) = MaxVals(0, X)
                MaxP(0) = X
            End If
            If MaxV(1) < MaxVals(1, X) Then
                MaxV(1) = MaxVals(1, X)
                MaxP(1) = X
            End If
        Next X
        If z = Target(0) Then
            CV(0, 0) = MaxV(0)
            CV(1, 0) = MaxV(1)
        End If
        If z = Target(1) Then
            CV(0, 1) = MaxV(0)
            CV(1, 1) = MaxV(1)
        End If
        MaxVals(0, MaxP(0)) = 0
        MaxVals(1, MaxP(1)) = 0
    Next z
    
    
    
    SS = GetTickCount
    ReDim PValMap(Len(StrainSeq(0)), PermNum)
    EE = GetTickCount
    TT = EE - SS
    X = X
    SS = GetTickCount
    Form1.ProgressBar1.Value = 90
    If X = X Then
        Dummy = MakePValMap(DN, Len(StrainSeq(0)), PermNum, Map(0, 0), PValMap(0, 0))
    Else
        Dim TopS As Double, BottomS As Double, DoneNum As Long
        For Y = DN To Len(StrainSeq(0)) - DN + 1
            TopS = 1000000: BottomS = 0
            DoneNum = 0
            For X = 0 To PermNum
                DoneNum = 0
                BottomS = 0
                For z = 0 To PermNum
                    
                    If Map(Y, z) > BottomS And Map(Y, z) < TopS Then
                        BottomS = Map(Y, z)
                        DoneNum = 1
                    End If
                Next z
                If DoneNum = 1 Then
                    DoneNum = 0
                    For z = 0 To PermNum
                        If Map(Y, z) = BottomS Then
                            PValMap(Y, X + DoneNum) = Map(Y, z)
                            DoneNum = DoneNum + 1
                        End If
                    Next z
                    
                    TopS = BottomS
                    
                    If DoneNum > 1 Then X = X + DoneNum - 1
                    LastX = X + 1
                Else
                    For z = LastX To PermNum
                        PValMap(Y, z) = 0
                    Next z
                    Exit For
                End If
            Next X
        Next Y
    End If
    EE = GetTickCount
    TT = EE - SS
    Form1.ProgressBar1.Value = 100
    X = X
End Sub
Public Sub DoMatCap()
    If CurMatrixFlag = 0 Then
        Form1.Label4.Caption = "Jakobsen's Compatibility Matrix"
    ElseIf CurMatrixFlag = 1 Then
        Form1.Label4.Caption = "Recombination Matrix"
    ElseIf CurMatrixFlag = 2 Then
        Form1.Label4.Caption = "Modulatrity Matrix"
    ElseIf CurMatrixFlag = 3 Then
        Form1.Label4.Caption = "Region Count Matrix"
    ElseIf CurMatrixFlag = 4 Then
        Form1.Label4.Caption = "Breakpoint Pair Matrix"
    ElseIf CurMatrixFlag = 5 Then
        Form1.Label4.Caption = "McVean's Linkage Disequilibrium Matrix"
    ElseIf CurMatrixFlag = 6 Then
        Form1.Label4.Caption = "Hudson and Kaplan's RMin Matrix"
    ElseIf CurMatrixFlag = 7 Then
        Form1.Label4.Caption = "Hudson and Kaplan's RMin/Distance Matrix"
    ElseIf CurMatrixFlag = 8 Then
        Form1.Label4.Caption = "MaxChi Breakpoint Matrix"
    ElseIf CurMatrixFlag = 11 Then
        Form1.Label4.Caption = "LARD Breakpoint Matrix"
    End If
End Sub
Public Sub ClearMatrix()
CurMatrixFlag = 255
Form1.Picture26.Picture = LoadPicture()
Form1.Picture18.Picture = LoadPicture()
Form1.Picture17.Picture = LoadPicture()
Form1.Label4.Caption = ""
                        
For z = 0 To 4
    Form1.Line1(z).Visible = False
    Form1.Label6(z) = ""
Next z
For z = 0 To 2
    Form1.Label7(z) = ""
Next z
Form1.VScroll5.Enabled = False
Form1.HScroll4.Enabled = False
End Sub

Public Sub DoHeatMaps()
    If X = X Then
        Y = -1
        For X = 128 To 255
            Y = Y + 1
            HeatMap(0, Y) = RGB(0, 0, X)
        Next X
        
        
        For X = 0 To 255
            
            HeatMap(0, X + 128) = RGB(0, X, 255)
        Next X
        
        Y = -1
        For X = 0 To 255 Step 2
            Y = Y + 1
            HeatMap(0, Y + 383) = RGB(0, 255, 255 - X)
        Next X
        Y = -1
        For X = 0 To 255 Step 2
            Y = Y + 1
            HeatMap(0, Y + 510) = RGB(X, 255, 0)
        Next X
        For X = 1 To 255
            HeatMap(0, X + 637) = RGB(255, 255 - X, 0)
        Next X
        Y = -1
        For X = 0 To 128
            Y = Y + 1
            HeatMap(0, Y + 893) = RGB(255 - X, 0, 0)
        Next X
        'X = X
    Else
        For X = 0 To 255
            HeatMap(0, X) = RGB(0, X, 255)
        Next X
        
        For X = 1 To 255
            HeatMap(0, X + 255) = RGB(0, 255, 255 - X)
        Next X
        
        For X = 1 To 255
            HeatMap(0, X + 510) = RGB(X, 255, 0)
        Next X
        For X = 1 To 255
            HeatMap(0, X + 765) = RGB(255, 255 - X, 0)
        Next X
    End If

'greyscale
    For X = 0 To 1020
        HeatMap(1, X) = RGB(255 - Int(X / 4), 255 - Int(X / 4), 255 - Int(X / 4))
    Next X
    For X = 0 To 1020
        HeatMap(2, X) = RGB(Int(X / 4), Int(X / 4), Int(X / 4))
    Next X
    'redscale
    For X = 0 To 510
        HeatMap(3, X) = RGB(Int(X / 2), 0, 0)
    Next X
    For X = 511 To 1020
        HeatMap(3, X) = RGB(255, 255 - (510 - X / 2), 255 - (510 - X / 2))
    Next X
    'greenscale
    For X = 0 To 510
        HeatMap(4, X) = RGB(0, Int(X / 2), 0)
    Next X
    For X = 511 To 1020
        HeatMap(4, X) = RGB(255 - (510 - X / 2), 255, 255 - (510 - X / 2))
    Next X
    'bluescale
    For X = 0 To 510
        HeatMap(5, X) = RGB(0, 0, Int(X / 2))
    Next X
    For X = 511 To 1020
        HeatMap(5, X) = RGB(255 - (510 - X / 2), 255 - (510 - X / 2), 255)
    Next X
    
    For X = 0 To 255
        HeatMap(6, X) = RGB(X, 0, 0)
    Next X
    For X = 0 To 510
        HeatMap(6, 256 + X) = RGB(255, X / 2, 0)
    Next X
    For X = 0 To 255
        HeatMap(6, 766 + X) = RGB(255, 255, X)
    Next X
End Sub
Public Sub DrawRecMatrix(PWFlag)
Dim Addj(1) As Double, AdjDst As Double, TParDist As Double, PairDist() As Double, PairValid As Double, PairDiff As Double, RegionMat() As Double, RSize As Long, SConvert As Double, ST As Long, EN As Long
Dim UseAll As Byte

    
SSS = GetTickCount

RSize = Len(StrainSeq(0))
If RSize > 2000 Then RSize = 2000
SConvert = RSize / Len(StrainSeq(0))
ReDim RegionMat(RSize + 1, RSize + 1)
    
If (PWFlag = 0 And DoneMatX(2) = 0) Or (PWFlag = 1 And DoneMatX(1) = 0) Then
    UseAll = 0
    Dim Size As Long, d As Long, P1 As Long, P2 As Long, Win As Long
        
    
    Dim PermutationX As Long
    PermutationX = 1000
    ENumb = SEventNumber
    
    Dim Excl() As Byte, Enu As Long, NC As Long
    ReDim Excl(SEventNumber), pCount(SEventNumber, AddNum)
        
    Dim BPV() As Double
    ReDim BPV(SEventNumber, AddNum)
    For X = 1 To SEventNumber
        For Y = 0 To AddNum
            BPV(X, Y) = LowestProb
        Next Y
    Next X
        
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            If XOverList(X, Y).ProgramFlag <= AddNum - 1 Then
                If BPV(SuperEventList(XOverList(X, Y).Eventnumber), XOverList(X, Y).ProgramFlag) > XOverList(X, Y).Probability And XOverList(X, Y).Probability > 0 Then
                    BPV(SuperEventList(XOverList(X, Y).Eventnumber), XOverList(X, Y).ProgramFlag) = XOverList(X, Y).Probability
                End If
            End If
        Next Y
    Next X
    
    Dim PValCon As Double
    For X = 1 To SEventNumber
        For Y = 0 To AddNum - 1
            If BPV(X, Y) = LowestProb And Confirm(X, Y) > 0 Then
                
                PValCon = ConfirmP(X, Y) / Confirm(X, Y)
                PValCon = 10 ^ (-PValCon)
                
                'ConfirmP(X, Y) = PValCon
                If BPV(X, Y) > PValCon Then
                    BPV(X, Y) = PValCon
                End If
            End If
        Next Y
    Next X
        
        
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            Enu = SuperEventList(XOverList(X, Y).Eventnumber)
                
            NC = 0
            If BPV(Enu, 0) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 1) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 3) < LowestProb Then
                    
                NC = NC + 1
                    
            ElseIf BPV(Enu, 4) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 2) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 5) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
                
            If XOverList(X, Y).Probability > 0 And XOverList(X, Y).Probability < LowestProb And NC > ConsensusProg And XOverList(X, Y).Accept <> 2 And XOverList(X, Y).ProgramFlag <= AddNum - 1 And XOverList(X, Y).MissIdentifyFlag <> 3 And XOverList(X, Y).MissIdentifyFlag <> 13 Then
                Excl(Enu) = 1
            
            End If
        Next Y
    Next X
        
        
    Dim BPos() As Long, APos() As Long, TypeSeq
    ReDim APos(Len(StrainSeq(0))), BPos(Len(StrainSeq(0)))
    
    TypeSeq = 0
    If TypeSeq > PermNextNo Then TypeSeq = 0
        
    For X = 1 To Len(StrainSeq(0))
        APos(X) = X - SeqSpaces(X, TypeSeq)
        BPos(X - SeqSpaces(X, TypeSeq)) = X
    Next X
    
    Win = 200
    If X = 12345 Then
        
        B = 0
        C = 0
        If BPCvalFlag = 0 Or X = X Then
            BPCvalFlag = 1
            ReDim BPCVal(1, 1)
            Dim DN As Long
            DN = 1
            Call RecombMapPermsB(DN, APos(), BPos(), Excl(), BPCVal(), Win, PermutationX, 3)
        End If
    End If
    
    zz = 0
            

    'Make exclusions
    Dim Age As Double, Tot(2) As Double, RS1 As Long, LS1 As Long, RS2 As Long, LS2 As Long, ID() As Double, NS(1) As Long, ParDist As Double
    xSeq1 = Seq1
    xSeq2 = Seq2
    xSeq3 = Seq3
    xrelx = RelX
    xrely = RelY
    xnjflag = NJFlag
    NJFlag = 0
    For G = 1 To ENumb
        
        If Excl(G) = 1 And (BestEvent(G, 0) > 0 Or BestEvent(G, 1) > 0) Then
            zz = zz + 1
            CNum = 0
            sen = X 'SuperEventlist(XOverList(BestEvent(X, 0), BestEvent(X, 1)).Eventnumber)
               
                
            d = XOverList(BestEvent(G, 0), BestEvent(G, 1)).Daughter
            P1 = XOverList(BestEvent(G, 0), BestEvent(G, 1)).MajorP
            P2 = XOverList(BestEvent(G, 0), BestEvent(G, 1)).MinorP
            ST = XOverList(BestEvent(G, 0), BestEvent(G, 1)).Beginning
            EN = XOverList(BestEvent(G, 0), BestEvent(G, 1)).Ending
            
            If XOverList(BestEvent(G, 0), BestEvent(G, 1)).Accept = 1 Or (UseAll = 1 And XOverList(BestEvent(G, 0), BestEvent(G, 1)).Accept <> 2) Then
                If PWFlag = 0 Or (PWFlag = 1 And (XOverList(BestEvent(G, 0), BestEvent(G, 1)).MissIdentifyFlag = 0 Or XOverList(BestEvent(G, 0), BestEvent(G, 1)).MissIdentifyFlag = 10)) Then
                
                    Seq1 = P1
                    Seq2 = P2
                    Seq3 = d
                    RelX = BestEvent(G, 0)
                    RelY = BestEvent(G, 1)
   
                    Call ModSeqNum(0)
                    
                    Call MakeTreeSeqs(ST, EN)
    
                    P1 = Seq1
                    P2 = Seq2
                    d = Seq3
                    'Work out the maximum identity of parental sequences at the time when the recombination event occured
                    
                    
                    ReDim ID(1, 2)
                    ID(0, 0) = FMat(d, P1): ID(0, 1) = FMat(d, P2): ID(0, 2) = FMat(P1, P2)
                    ID(1, 0) = SMat(d, P1): ID(1, 1) = SMat(d, P2): ID(1, 2) = SMat(P1, P2)
                    'Work out the distance modifyer for each region
                    Tot(0) = 0: Tot(1) = 0: Tot(2) = 0
                    '@'@'@'@'@
                    If X = X Then
                        TParDist = MakeDistModPar(Nextno, UBound(FMat, 1), UBound(PermDiffs, 1), Tot(0), FMat(0, 0), SMat(0, 0), PermDiffs(0, 0), PermValid(0, 0))
                        
                    Else
                        For X = 0 To Nextno
                            For Y = X + 1 To Nextno
                                If FMat(X, Y) < 3 Then
                                    Tot(0) = Tot(0) + FMat(X, Y)
                                    Tot(1) = Tot(1) + SMat(X, Y)
                                    TParDist = 1 - PermDiffs(X, Y) / PermValid(X, Y)
                                    If TParDist > 0.25 Then
                                        TParDist = (4# * TParDist - 1#) / 3#
                                        TParDist = Log(TParDist)
                                        TParDist = -0.75 * TParDist
                                       
                                    Else
                                        TParDist = 1
                                    End If
                                    Tot(2) = Tot(2) + TParDist
                                End If
                            Next Y
                        Next X
                    End If
                    If Tot(0) > 0 And Tot(2) > 0 Then '42.6091706,18.498428
                        Addj(0) = Tot(0) / Tot(2)
                    Else
                        Addj(0) = 1
                    End If
                    If Tot(1) > 0 And Tot(2) > 0 Then '42.6091706,18.498428
                        Addj(1) = Tot(1) / Tot(2)
                    Else
                        Addj(1) = 1
                    End If
                    'Modify the recombinant region dists
                    For X = 0 To 2
                        ID(0, X) = ID(0, X) / Addj(0)
                        ID(1, X) = ID(1, X) / Addj(1)
                        X = X
                    Next X
                    
                    ParDist = 0: Age = 1000
                    For X = 0 To 1
                        For Y = 0 To 2
                            If Age > ID(X, Y) Then Age = ID(X, Y)
                            If ParDist < ID(X, Y) Then ParDist = ID(X, Y)
                            XX = FMat(d, P1)
                        Next Y
                    Next X
                    ParDist = ParDist - Age
                    If ParDist < 0 Then ParDist = 0
                    If PWFlag = 1 Then
                        Dim PUse() As Long
                        ReDim PUse(Len(StrainSeq(0)))
                        If ST < EN Then
                            For z = 1 To ST - 1
                                PUse(z) = P2
                            Next z
                            For z = ST To EN
                                PUse(z) = P1
                            Next z
                            For z = EN + 1 To Len(StrainSeq(0))
                                PUse(z) = P2
                            Next z
                        Else
                            For z = 1 To EN
                                PUse(z) = P1
                            Next z
                            For z = EN + 1 To ST - 1
                                PUse(z) = P2
                            Next z
                            For z = ST To Len(StrainSeq(0))
                                PUse(z) = P1
                            Next z
                        End If
                        
                        pwin = 25
                        ReDim PairDist(Len(StrainSeq(0)))
                        PairValid = 0: PairDiff = 0
                        For z = 1 - pwin To pwin + 1
                            
                            If z < 1 Then
                                X = z + Len(StrainSeq(0))
                            ElseIf z > Len(StrainSeq(0)) Then
                                X = z - Len(StrainSeq(0))
                            Else
                                X = z
                            End If
                            If SeqNum(X, d) <> 46 And SeqNum(X, PUse(X)) <> 46 Then
                                PairValid = PairValid + 1
                                If SeqNum(X, d) <> SeqNum(X, PUse(X)) Then
                                    PairDiff = PairDiff + 1
                                End If
                                
                            End If
                        Next z
                        If PairValid > 0 Then
                            PairDist(1) = PairDiff / PairValid
                        Else
                            PairDist(1) = 0
                        End If
                        For Y = 2 To Len(StrainSeq(0))
                            z = Y - pwin
                            If z < 1 Then
                                X = z + Len(StrainSeq(0))
                            ElseIf z > Len(StrainSeq(0)) Then
                                X = z - Len(StrainSeq(0))
                            Else
                                X = z
                            End If
                            If SeqNum(X, d) <> 46 And SeqNum(X, PUse(X)) <> 46 Then
                                
                                PairValid = PairValid - 1
                                If SeqNum(X, d) <> SeqNum(X, PUse(X)) Then
                                    PairDiff = PairDiff - 1
                                End If
                            
                            End If
                            
                            z = Y + pwin
                            If z < 1 Then
                                X = z + Len(StrainSeq(0))
                            ElseIf z > Len(StrainSeq(0)) Then
                                X = z - Len(StrainSeq(0))
                            Else
                                X = z
                            End If
                            
                            If SeqNum(X, d) <> 46 And SeqNum(X, PUse(X)) <> 46 Then
                                
                                PairValid = PairValid + 1
                                If SeqNum(X, d) <> SeqNum(X, PUse(X)) Then
                                    PairDiff = PairDiff + 1
                                End If
                            
                            End If
                            If PairValid > 0 Then
                                PairDist(Y) = PairDiff / PairValid '1194,1195,1196,1197,1198 =3/12; 1199,1200 = 4/12
                                
                            Else
                                PairDist(Y) = 0
                            End If
                        Next Y
                        XX = UBound(PermValid, 1)
                        If P1 <= UBound(PermValid, 1) And P2 <= UBound(PermValid, 1) Then
                            If PermValid(P1, P2) + SubValid(P1, P2) > 0 Then
                                TParDist = 1 - PermDiffs(P1, P2) / PermValid(P1, P2)
                                If TParDist > 0.25 Then
                                    TParDist = (4# * TParDist - 1#) / 3#
                                    TParDist = Log(TParDist)
                                    TParDist = -0.75 * TParDist
                                   
                                Else
                                    TParDist = 1
                                End If
                            Else
                                TParDist = 0
                            End If
                        Else
                            TParDist = 0
                        End If
                        If TParDist > 0 Then
                            AdjDst = ParDist / TParDist
                            If AdjDst > 1 Then
                                AdjDst = 1
                               
                            End If
                            For X = 1 To Len(StrainSeq(0))
                                
                                PairDist(X) = PairDist(X) * AdjDst
                                
                                If PairDist(X) > 0.5 Then PairDist(X) = 0.5
                            Next X
                        Else
                            AdjDst = 0
                        End If
                        LS1 = -1: LS2 = -1
                            If ST < EN Then
                                For A = 1 To ST - 1
                                    RS1 = CInt(A * SConvert)
                                    If LS1 <> RS1 Then
                                        LS1 = RS1
                                    
                                        For B = ST To EN
                                            RS2 = CInt(B * SConvert)
                                            If LS2 <> RS2 Then
                                                LS2 = RS2
                                                
                                                If PairDist(A) > PairDist(B) Then
                                                    mdst = PairDist(B)
                                                Else
                                                    mdst = PairDist(A)
                                                End If
                                                If mdst > RegionMat(RS2, RS1) Then
                                                    RegionMat(RS2, RS1) = mdst
                                                    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                End If
                                                'If (PairDist(A) + PairDist(B)) / 2 > RegionMat(RS2, RS1) Then
                                                '    RegionMat(RS2, RS1) = (PairDist(A) + PairDist(B)) / 2
                                                '    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                'End If
                                            End If
                                        Next B
                                    End If
                                Next A
                                For A = EN + 1 To Len(StrainSeq(0))
                                    RS1 = CInt(A * SConvert)
                                    If LS1 <> RS1 Then
                                        LS1 = RS1
                                    
                                        For B = ST To EN
                                            RS2 = CInt(B * SConvert)
                                            If LS2 <> RS2 Then
                                                LS2 = RS2
                                                If PairDist(A) > PairDist(B) Then
                                                    mdst = PairDist(B)
                                                Else
                                                    mdst = PairDist(A)
                                                End If
                                                If mdst > RegionMat(RS2, RS1) Then
                                                    RegionMat(RS2, RS1) = mdst
                                                    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                End If
                                                'If (PairDist(A) + PairDist(B)) / 2 > RegionMat(RS2, RS1) Then
                                                '    RegionMat(RS2, RS1) = (PairDist(A) + PairDist(B)) / 2
                                                '    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                'End If
                                            End If
                                        Next B
                                    End If
                                Next A
                            Else
                                For A = EN + 1 To ST - 1
                                    RS1 = CInt(A * SConvert)
                                    If LS1 <> RS1 Then
                                        LS1 = RS1
                                    
                                        For B = 1 To EN
                                            RS2 = CInt(B * SConvert)
                                            If LS2 <> RS2 Then
                                                LS2 = RS2
                                                If PairDist(A) > PairDist(B) Then
                                                    mdst = PairDist(B)
                                                Else
                                                    mdst = PairDist(A)
                                                End If
                                                If mdst > RegionMat(RS2, RS1) Then
                                                    RegionMat(RS2, RS1) = mdst
                                                    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                End If
                                                'If (PairDist(A) + PairDist(B)) / 2 > RegionMat(RS2, RS1) Then
                                                '    RegionMat(RS2, RS1) = (PairDist(A) + PairDist(B)) / 2
                                                '    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                'End If
                                            End If
                                        Next B
                                        For B = ST To Len(StrainSeq(0))
                                            RS2 = CInt(B * SConvert)
                                            If LS2 <> RS2 Then
                                                LS2 = RS2
                                                If PairDist(A) > PairDist(B) Then
                                                    mdst = PairDist(B)
                                                Else
                                                    mdst = PairDist(A)
                                                End If
                                                If mdst > RegionMat(RS2, RS1) Then
                                                    RegionMat(RS2, RS1) = mdst
                                                    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                End If
                                                'If (PairDist(A) + PairDist(B)) / 2 > RegionMat(RS2, RS1) Then
                                                '    RegionMat(RS2, RS1) = (PairDist(A) + PairDist(B)) / 2
                                                '    RegionMat(RS1, RS2) = RegionMat(RS2, RS1)
                                                'End If
                                            End If
                                        Next B
                                    End If
                                Next A
                                
                            End If
                    Else
                    
                    
                        If ParDist > 0 Then
                            '@'@'@'@'@'@'@'@'@
                            If X = X Then
                                Dummy = AddToRegionMat(Len(StrainSeq(0)), ST, EN, RSize, SConvert, ParDist, RegionMat(0, 0))
                            Else
                                LS1 = -1: LS2 = -1
                                If ST < EN Then
                                    
                                    For A = 1 To ST - 1
                                        RS1 = CInt(A * SConvert)
                                        If LS1 <> RS1 Then
                                            LS1 = RS1
                                            
                                            For B = ST To EN
                                                RS2 = CInt(B * SConvert)
                                                If LS2 <> RS2 Then
                                                    LS2 = RS2
                                                    
                                                    If ParDist > RegionMat(RS2, RS1) Then
                                                        RegionMat(RS2, RS1) = ParDist
                                                        RegionMat(RS1, RS2) = ParDist
                                                    End If
                                                End If
                                            Next B
                                        End If
                                    Next A
                                   
                                    For A = EN + 1 To Len(StrainSeq(0))
                                        RS1 = CInt(A * SConvert)
                                        If LS1 <> RS1 Then
                                            LS1 = RS1
                                        
                                            For B = ST To EN
                                                RS2 = CInt(B * SConvert)
                                                If LS2 <> RS2 Then
                                                    LS2 = RS2
                                                    
                                                    If ParDist > RegionMat(RS2, RS1) Then
                                                        RegionMat(RS2, RS1) = ParDist
                                                        RegionMat(RS1, RS2) = ParDist
                                                    End If
                                                End If
                                            Next B
                                        End If
                                    Next A
                                Else
                                    
                                    For A = EN + 1 To ST - 1
                                        RS1 = CInt(A * SConvert)
                                        If LS1 <> RS1 Then
                                            LS1 = RS1
                                        
                                            For B = 1 To EN
                                                RS2 = CInt(B * SConvert)
                                                If LS2 <> RS2 Then
                                                    LS2 = RS2
                                                    
                                                    If ParDist > RegionMat(RS2, RS1) Then
                                                        RegionMat(RS2, RS1) = ParDist
                                                        RegionMat(RS1, RS2) = ParDist
                                                    End If
                                                End If
                                            Next B
                                            For B = ST To Len(StrainSeq(0))
                                                RS2 = CInt(B * SConvert)
                                                If LS2 <> RS2 Then
                                                    LS2 = RS2
                                                    If ParDist > RegionMat(RS2, RS1) Then
                                                        RegionMat(RS2, RS1) = ParDist
                                                        RegionMat(RS1, RS2) = ParDist
                                                    End If
                                                End If
                                            Next B
                                        End If
                                    Next A
                                    
                                End If
                            End If
                        End If
                    End If
                    Call UnModSeqNum(0)
                End If
            End If
    
            
        End If
        
        Form1.SSPanel1.Caption = Trim(Str(G)) + " of " + Trim(Str(SEventNumber)) + " events mapped"
        Form1.SSPanel1.Refresh
        Form1.ProgressBar1 = (G / SEventNumber) * 70
        Form1.Refresh
    Next G
    
    
    
    Seq1 = xSeq1
    Seq2 = xSeq2
    Seq3 = xSeq3
    RelX = xrelx
    RelY = xrely
    NJFlag = xnjflag
Else
    SS = GetTickCount
    If X = 12345 Then
        If PWFlag = 1 Then
            For X = 0 To RSize + 1
                For Y = 0 To RSize + 1
                    RegionMat(X, Y) = MatrixM(X, Y)
                Next Y
            Next X
        Else
            For X = 0 To RSize + 1
                For Y = 0 To RSize + 1
                    RegionMat(X, Y) = MatrixR(X, Y)
                Next Y
            Next X
        End If
    End If
    EE = GetTickCount
    TT = EE - SS
    X = X
End If

Dim MaxN As Double, MinN As Double






If (PWFlag = 0 And DoneMatX(2) = 0) Or (PWFlag = 1 And DoneMatX(1) = 0) Then
    MaxN = FindMaxN(RSize, RegionMat(0, 0))
    If PWFlag = 0 Then
        MatBound(1) = MaxN
        ReDim MatrixR(RSize + 1, RSize + 1)
        If X = X Then
            UB1 = UBound(RegionMat, 1): UB2 = UBound(RegionMat, 2)
            UB3 = UBound(MatrixR, 1): UB4 = UBound(MatrixR, 2)
            
            Dummy = CopyDoubleArray(RSize + 1, RSize + 1, UB1, UB2, UB3, UB4, RegionMat(0, 0), MatrixR(0, 0))
        Else
        
            For X = 0 To RSize + 1
                For Y = 0 To RSize + 1
                    MatrixR(X, Y) = RegionMat(X, Y)
                Next Y
            Next X
        End If
    Else
        MatBound(2) = MaxN
        
        
        ReDim MatrixM(RSize + 1, RSize + 1)
        If X = X Then
            UB1 = UBound(RegionMat, 1): UB2 = UBound(RegionMat, 2)
            UB3 = UBound(MatrixM, 1): UB4 = UBound(MatrixM, 2)
            
            Dummy = CopyDoubleArray(RSize + 1, RSize + 1, UB1, UB2, UB3, UB4, RegionMat(0, 0), MatrixM(0, 0))
        Else
        
            For X = 0 To RSize + 1
                For Y = 0 To RSize + 1
                    MatrixM(X, Y) = RegionMat(X, Y)
                Next Y
            Next X
        End If
    End If

    If MaxN = 0 Then
        
        Form1.ProgressBar1 = 0
        Form1.SSPanel1 = ""
        Exit Sub
    End If

    Form1.SSPanel1.Caption = "Drawing matrix"
Else
    If PWFlag = 0 Then
        MaxN = MatBound(1)
    Else
        MaxN = MatBound(2)
    End If
End If




    Form1.Picture26.ScaleMode = 3
    Form1.Picture26.AutoRedraw = True
    
    Dim XAddj As Double, PosS(1) As Long, PosE(1) As Long, DistD As Long
    If PWFlag = 0 Then
        DistD = RSize / MatZoom(1)
        PosS(0) = MatCoord(1, 0)
        PosE(0) = PosS(0) + DistD
        PosS(1) = MatCoord(1, 1)
        PosE(1) = PosS(1) + DistD
        If PosE(1) > (UBound(MatrixR, 1) - 1) - 1 Then PosE(1) = (UBound(MatrixR, 1) - 1) - 1
        If PosE(0) > (UBound(MatrixR, 1) - 1) - 1 Then PosE(0) = (UBound(MatrixR, 1) - 1) - 1
       
        If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
        If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
    Else
        DistD = RSize / MatZoom(2)
        PosS(0) = MatCoord(2, 0)
        PosE(0) = PosS(0) + DistD
        PosS(1) = MatCoord(2, 1)
        PosE(1) = PosS(1) + DistD
        If PosE(1) > (UBound(MatrixM, 1) - 1) - 1 Then PosE(1) = (UBound(MatrixM, 1) - 1) - 1
        If PosE(0) > (UBound(MatrixM, 1) - 1) - 1 Then PosE(0) = (UBound(MatrixM, 1) - 1) - 1
       
        If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
        If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
        
    End If
    Form1.Picture26.Picture = LoadPicture()
    XAddj = (Form1.Picture26.ScaleHeight) / DistD
    
    'For B = 1 To 20
        Pict = Form1.Picture26.hDC
        
        If X = 1234 Then
            Dummy = DrawMats(Pict, MaxN, CurScale, XAddj, UBound(RegionMat, 1), PosS(0), PosE(0), RegionMat(0, 0), HeatMap(0, 0))
            
        Else
            If PWFlag = 0 Then
                Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixR(), HeatMap(), CurScale, MaxN)
            Else
                Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixM(), HeatMap(), CurScale, MaxN)
            End If
        End If
        
    'Next B
    Form1.Picture26.Refresh
    Form1.ProgressBar1 = 0
    Form1.SSPanel1 = ""
    If PWFlag = 0 Then
        DoneMatX(2) = 1
        If DontDoKey = 0 Then
            Call DoKey(0, MaxN, MinN, 1, "Average parental genetic distance", CurScale) ', "relatedness")
        End If
    Else
        DoneMatX(1) = 1
        If DontDoKey = 0 Then
            
            Call DoKey(0, MaxN, MinN, 2, "Parental genetic distance", CurScale) ', "relatedness")
        End If
    End If
    
    Call UnModNextno
    EE = GetTickCount
    TT = EE - SSS
    X = X
    '121 seconds
    '95 seconds
    '50.104 seconds
    '31.8 seconds
    '30.844
    '28.547
    '26.531
    '22.625
    '20.812
    '19.063
    
End Sub
Public Sub DrawmatsVB(Min, PosE() As Long, PosS() As Long, SX, SY, XAddj, RegionMat() As Double, HeatMap() As Long, CurScale, MaxN)
    Dim StS As Long, RM As Double, PntAPI As POINTAPI, YP As Long, z As Long, X As Long, Y As Long, Pict As Long, StSX
    Form1.Picture26.ScaleMode = 3
    
    'XAddj = (Form1.Picture26.ScaleHeight) / DistD
    Form1.Picture26.AutoRedraw = True
    SS = GetTickCount
    Pict = Form1.Picture26.hDC
    StS = (1 / XAddj) - 1
    If StS < 1 Then StS = 1
    StSX = StS
    If StSX < 1 Then StSX = 1
    Dim MR
    MR = MaxN - Min
    Dim XD As Double, ZD As Double
    If XAddj <= 1 Then
        For Y = SY To PosE(1) Step StS
                    
                For X = SX To PosE(0) Step StSX
                        
                            
                        z = X + StSX
                        RM = RegionMat(X, Y)
                        Do While z <= PosE(0)
                            If RegionMat(z, Y) <> RM Then Exit Do
                            z = z + StSX
                        Loop
                        z = z - StSX
                        
                        
                        If z = X Then 'and Then ' And X = 12345 Then
                            SetPixelV Pict, Int((X - PosS(0)) * XAddj), Int((Y - PosS(1)) * XAddj), HeatMap(CurScale, CInt(((RegionMat(X, Y) - Min) / MR) * 1020))
                        Else
                            YP = Int((Y - PosS(1)) * XAddj)
                            
                            Form1.Picture26.ForeColor = HeatMap(CurScale, CInt(((RegionMat(X, Y) - Min) / MR) * 1020))
                            Dummy = MoveToEx(Pict, Int((X - PosS(0)) * XAddj), YP, PntAPI)
                            Dummy = LineTo(Pict, Int((z - PosS(0)) * XAddj) + 1, YP)
                            X = z
                        End If
                        
                Next X
        Next Y
    Else
        For Y = SY To PosE(1) Step StS
                    
                For X = SX To PosE(0) Step StSX
                        
                            
                        z = X + StSX
                        RM = RegionMat(X, Y)
                        Do While z <= PosE(0)
                            If RegionMat(z, Y) <> RM Then Exit Do
                            z = z + StSX
                        Loop
                        z = z - StSX
                        ZD = z + 0.5
                        XD = X - 0.5
                        
                        'If Z = X Then 'and Then ' And X = 12345 Then
                        '    SetPixelV Pict, Int((X - PosS(0)) * XAddj), Int((Y - PosS(1)) * XAddj), HeatMap(CurScale, CInt(((RegionMat(X, Y) - Min) / MR) * 1020))
                        'Else
                            band = (CLng(XAddj) + 1) / 2
                           
                            
                            YP = Int((Y - PosS(1)) * XAddj)
                            syp = YP - band
                            eyp = YP + band
                            For YP = syp To eyp
                                Form1.Picture26.ForeColor = HeatMap(CurScale, CInt(((RegionMat(X, Y) - Min) / MR) * 1020))
                                Dummy = MoveToEx(Pict, Int((XD - PosS(0)) * XAddj), YP, PntAPI)
                                Dummy = LineTo(Pict, Int((ZD - PosS(0)) * XAddj) + 1, YP)
                            Next YP
                            X = z
                            
                        'End If
                        
                Next X
        Next Y
    End If
    EE = GetTickCount
    TT = EE - SSS
    X = X
    '121.500
    '95.828
End Sub
Public Sub DrawBPMatrix()
Dim RSize As Long, SConvert As Double

RSize = Len(StrainSeq(0))
If RSize > 1000 Then RSize = 1000
SConvert = RSize / Len(StrainSeq(0))

Dim UseAll As Byte
UseAll = 1
If DoneMatX(4) = 0 Then
    ReDim MatrixBP(RSize + 1, RSize + 1)
    Dim Size As Long, d As Long, P1 As Long, P2 As Long, Win As Long
        
    
    Dim PermutationX As Long
    PermutationX = 1000
    ENumb = SEventNumber
    
    Dim Excl() As Byte, Enu As Long, NC As Long
    ReDim Excl(SEventNumber), pCount(SEventNumber, AddNum)
        
    Dim BPV() As Double
    ReDim BPV(SEventNumber, AddNum)
    For X = 1 To SEventNumber
        For Y = 0 To AddNum
            BPV(X, Y) = LowestProb
        Next Y
    Next X
        
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            If XOverList(X, Y).ProgramFlag <= AddNum - 1 Then
                If BPV(SuperEventList(XOverList(X, Y).Eventnumber), XOverList(X, Y).ProgramFlag) > XOverList(X, Y).Probability And XOverList(X, Y).Probability > 0 Then
                    BPV(SuperEventList(XOverList(X, Y).Eventnumber), XOverList(X, Y).ProgramFlag) = XOverList(X, Y).Probability
                End If
            End If
        Next Y
    Next X
    
    Dim PValCon As Double
    For X = 1 To SEventNumber
        For Y = 0 To AddNum - 1
            If BPV(X, Y) = LowestProb And Confirm(X, Y) > 0 Then
                
                PValCon = ConfirmP(X, Y) / Confirm(X, Y)
                PValCon = 10 ^ (-PValCon)
                
                'ConfirmP(X, Y) = PValCon
                If BPV(X, Y) > PValCon Then
                    BPV(X, Y) = PValCon
                End If
            End If
        Next Y
    Next X
        
        
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            Enu = SuperEventList(XOverList(X, Y).Eventnumber)
           ' If Enu = 2 Then
           '     X = X
           ' End If
            NC = 0
            If BPV(Enu, 0) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 1) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 3) < LowestProb Then
                    
                NC = NC + 1
                    
            ElseIf BPV(Enu, 4) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 2) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 5) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
                
            If XOverList(X, Y).Probability > 0 And XOverList(X, Y).Probability < LowestProb And NC > ConsensusProg And XOverList(X, Y).Accept <> 2 And XOverList(X, Y).ProgramFlag <= AddNum - 1 And ((XOverList(X, Y).MissIdentifyFlag <> 3 And XOverList(X, Y).MissIdentifyFlag <> 13) Or (XOverList(X, Y).Accept = 1)) Then
                Excl(Enu) = 1
            
            End If
            X = X
        Next Y
    Next X
        
        
    Dim BPos() As Long, APos() As Long, TypeSeq
    ReDim APos(Len(StrainSeq(0))), BPos(Len(StrainSeq(0)))
    
    TypeSeq = 0
    If TypeSeq > PermNextNo Then TypeSeq = 0
        
    For X = 1 To Len(StrainSeq(0))
        APos(X) = X - SeqSpaces(X, TypeSeq)
        BPos(X - SeqSpaces(X, TypeSeq)) = X
    Next X
        
    If X = 12345 Then
        Win = 200
        B = 0
        C = 0
        If BPCvalFlag = 0 Or X = X Then
            BPCvalFlag = 1
            ReDim BPCVal(1, 1)
            Dim DN As Long
            DN = 1
            Call RecombMapPermsB(DN, APos(), BPos(), Excl(), BPCVal(), Win, PermutationX, 3)
        End If
    End If
    
    zz = 0
    Win = MatWinSize / 2
    'Make exclusions
    Dim RS1 As Long, LS1 As Long, RS2 As Long, LS2 As Long
    For X = 1 To ENumb
        
        If Excl(X) = 1 And (BestEvent(X, 0) > 0 Or BestEvent(X, 1) > 0) And XOverList(BestEvent(X, 0), BestEvent(X, 1)).SBPFlag = 0 Then
            zz = zz + 1
            CNum = 0
            sen = X 'SuperEventlist(XOverList(BestEvent(X, 0), BestEvent(X, 1)).Eventnumber)
               
                
            d = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Daughter
            P1 = XOverList(BestEvent(X, 0), BestEvent(X, 1)).MajorP
            P2 = XOverList(BestEvent(X, 0), BestEvent(X, 1)).MinorP
            ST = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Beginning
            EN = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Ending
            LS1 = -1: LS2 = -1
            For A = EN - Win + 1 To EN + Win
                    
                    RS1 = CInt(A * SConvert)
                    
                    If RS1 < 1 Then RS1 = RSize + RS1
                    If RS1 > RSize Then RS1 = RS1 - RSize
                    If LS1 <> RS1 Then
                        
                        LS1 = RS1
                    
                        For B = ST - Win + 1 To ST + Win
                            If Sqr((A - EN) ^ 2 + (B - ST) ^ 2) <= Win Then
                                RS2 = CInt(B * SConvert)
                                If RS2 < 1 Then RS2 = RSize + RS2
                                If RS2 > RSize Then RS2 = RS2 - RSize
                                If LS2 <> RS2 Then
                                    
                                    LS2 = RS2
                                    MatrixBP(RS2, RS1) = MatrixBP(RS2, RS1) + 1
                                    MatrixBP(RS1, RS2) = MatrixBP(RS1, RS2) + 1
                                End If
                            End If
                        Next B
                    End If
            Next A
                
            
    
            
        End If
        
        Form1.SSPanel1.Caption = Trim(Str(X)) + " of " + Trim(Str(SEventNumber)) + " breakpoints mapped"
        Form1.ProgressBar1 = (X / SEventNumber) * 100
    Next X
    Form1.ProgressBar1 = 100
    DoneMatX(4) = 1
    
     
        Form1.Picture26.ScaleMode = 3
        'Form1.Picture26.Picture = LoadPicture()
        Form1.Picture26.AutoRedraw = True
       ' Form1.Picture26.ScaleHeight = Form1.Picture26.ScaleWidth
        
        Dim XAddj As Double
        
        
        Dim MaxN As Double
            MaxN = FindMaxN(RSize, MatrixBP(0, 0))
            'MaxN = MaxN + 1
            MatBound(4) = MaxN
        Form1.SSPanel1.Caption = "Drawing matrix"
    Else
        MaxN = MatBound(4)
    End If
    Form1.Picture26.Picture = LoadPicture()
    
    
    Dim PosS(1) As Long, PosE(1) As Long, DistD As Long
        XAddj = (Form1.Picture26.ScaleHeight) / RSize
        DistD = RSize / MatZoom(4)
        PosS(0) = MatCoord(4, 0)
        PosE(0) = PosS(0) + DistD
        PosS(1) = MatCoord(4, 1)
        PosE(1) = PosS(1) + DistD
        If PosE(1) > (UBound(MatrixBP, 1) - 1) - 1 Then PosE(1) = (UBound(MatrixBP, 1) - 1) - 1
        If PosE(0) > (UBound(MatrixBP, 1) - 1) - 1 Then PosE(0) = (UBound(MatrixBP, 1) - 1) - 1
       
        If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
        If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
        
        If MaxN > 4 Then
            MaxN = CLng((MaxN / 4) + 0.5) * 4
        ElseIf MaxN >= 4 Then
            MaxN = CLng((MaxN / 2) + 0.5) * 2
        Else
            MaxN = 2
        End If
        Form1.Picture26.Picture = LoadPicture()
    XAddj = (Form1.Picture26.ScaleHeight) / DistD
    If X = 12345 Then
        Pict = Form1.Picture26.hDC
        Dummy = MakeHeatPlot(MaxN, UBound(HeatMap, 1), CurScale, Pict, RSize, XAddj, HeatMap(0, 0), MatrixBP(0, 0))
        
    Else
        If X = X Then
            Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixBP(), HeatMap(), CurScale, MaxN)
        Else
    
            RangeN = MaxN - MinN
        'For B = 1 To 20
            SS = GetTickCount
            Pict = Form1.Picture26.hDC
            For X = 1 To RSize 'Step 20
                For Y = X To RSize 'Step 20
                    'Form1.Picture26.PSet (CInt(X * XAddj), CInt(Y * XAddj)), heatmap(curscale,CInt((matrixbp(X, Y) / MaxN) * 1020))
                    'Form1.Picture26.PSet (CInt(Y * XAddj), CInt(X * XAddj)), heatmap(curscale,CInt((matrixbp(X, Y) / MaxN) * 1020))
                    SetPixelV Pict, CInt(X * XAddj), CInt(Y * XAddj), HeatMap(CurScale, CInt((MatrixBP(X, Y) / MaxN) * 1020))
                    SetPixelV Pict, CInt(Y * XAddj), CInt(X * XAddj), HeatMap(CurScale, CInt((MatrixBP(X, Y) / MaxN) * 1020))
                Next Y
                Form1.Picture26.Refresh
                Form1.ProgressBar1 = 20 + (X / RSize) * 80
            Next X
        End If
    End If
    If DontDoKey = 0 Then
        Call DoKey(0, MaxN, MinN, 4, "Number of breakpoints", CurScale)
    End If
    Form1.Picture26.Refresh
        EE = GetTickCount
        TT = EE - SS
        X = X
    'Next B
    Form1.ProgressBar1 = 0
    Form1.SSPanel1 = ""
End Sub
Public Sub DoKey(DiscF, MaxN, MinN, TF, Cap1 As String, CScale)

Dim PH As Long, LongCap As Long, ML As Long

XX = X
For X = 0 To 4
    XX = Form1.Line1(X).Visible
    Form1.Line1(X).Visible = False
Next X

'Form1.Picture17.Top = Form1.Picture29.Top + 10 * Screen.TwipsPerPixelY
'Exit Sub


Dim tMaxN As Long
tMaxN = MaxN

'MaxN = 1
If MaxN >= 3 Or TF = 1 Or TF = 2 Or DiscF = 1 Then
    Form1.Label6(0).Caption = CLng((MaxN) * 1000) / 1000
    Form1.Label6(1).Caption = CLng((MaxN - (MaxN - MinN) / 4) * 1000) / 1000
    Form1.Label6(2).Caption = CLng((MaxN - (MaxN - MinN) / 2) * 1000) / 1000
    Form1.Label6(3).Caption = CLng((MinN + (MaxN - MinN) / 4) * 1000) / 1000
    Form1.Label6(4).Caption = CLng(MinN * 1000) / 1000
    For X = 0 To 4
        Form1.Line1(X).Visible = True
    Next X

ElseIf MaxN = 2 Then
    Form1.Label6(0).Caption = CLng((MaxN) * 1000) / 1000
    Form1.Label6(1).Caption = "" 'CLng((MaxN - (MaxN - MinN) / 4) * 1000) / 1000
    Form1.Label6(2).Caption = CLng((MaxN - (MaxN - MinN) / 2) * 1000) / 1000
    Form1.Label6(3).Caption = "" 'CLng((MinN + (MaxN - MinN) / 4) * 1000) / 1000
    Form1.Label6(4).Caption = CLng(MinN * 1000) / 1000
    
    Form1.Line1(0).Visible = True
    Form1.Line1(2).Visible = True
    Form1.Line1(4).Visible = True
    
ElseIf MaxN = 1 Then
    Form1.Label6(0).Caption = ""
    Form1.Label6(1).Caption = "1"
    Form1.Label6(2).Caption = ""
    Form1.Label6(3).Caption = "0"
    Form1.Label6(4).Caption = ""
    
    Form1.Line1(1).Visible = True
    Form1.Line1(3).Visible = True
End If
ML = 0: LongCap = 0
For X = 0 To 4
    If Len(Form1.Label6(X).Caption) > ML Then
        ML = Len(Form1.Label6(X).Caption)
        LongCap = X
        
    End If
Next X
'
If Form1.Picture17.Top <> Form1.Picture29.Top + 10 * Screen.TwipsPerPixelY Then
    Diff = Form1.Picture17.Top - (Form1.Picture29.Top + 10 * Screen.TwipsPerPixelY)
    Form1.Picture17.Top = Form1.Picture17.Top - Diff
    For X = 0 To 4
        Form1.Line1(X).Y1 = Form1.Line1(X).Y1 - Diff
        Form1.Line1(X).Y2 = Form1.Line1(X).Y2 - Diff
        Form1.Label6(X).Top = Form1.Label6(X).Top - Diff
    Next X
    
End If


Form1.Picture18.AutoRedraw = 1
Form1.Picture18.Picture = LoadPicture()

Form1.Picture18.ScaleMode = 3
Form1.Picture18.CurrentX = 11 '(Form1.Picture18.Width / Screen.TwipsPerPixelY) - 5
Form1.Picture18.FontSize = 5
TW = ((Form1.Picture18.TextWidth(Cap1) * 1.3) / Screen.TwipsPerPixelX) * Screen.TwipsPerPixelY
Form1.Picture18.Height = TW * Screen.TwipsPerPixelY
'If TW > Form1.Picture18.Height / Screen.TwipsPerPixelY Then
'    Form1.Picture18.Height = TW * Screen.TwipsPerPixelY
'
'End If
Form1.Picture18.Top = Form1.Picture17.Top
Form1.Picture18.Top = Form1.Picture18.Top - (Form1.Picture18.Height - Form1.Picture17.Height) / 2
Form1.Picture18.CurrentY = ((Form1.Picture18.Height / Screen.TwipsPerPixelY) - TW) / 2
Form1.Picture18.left = Form1.Label6(LongCap).left + Form1.Label6(LongCap).Width + 50 '+ ((Form1.SSPanel6(2).Width - 50)-(Form1.Label6(LongCap).left + Form1.Label6(LongCap).Width + 50)/2 '((Form1.SSPanel6(2).Width - 50 - Form1.Picture18.Width) - Form1.Label6(LongCap).left + Form1.Label6(LongCap).Width) / 2
If Form1.Picture18.Top < Form1.Picture29.Top Then
    Diff = Form1.Picture29.Top - Form1.Picture18.Top
    Form1.Picture17.Top = Form1.Picture17.Top + Diff
    Form1.Picture18.Top = Form1.Picture18.Top + Diff
    For X = 0 To 4
        Form1.Line1(X).Y1 = Form1.Line1(X).Y1 + Diff
        Form1.Line1(X).Y2 = Form1.Line1(X).Y2 + Diff
        Form1.Label6(X).Top = Form1.Label6(X).Top + Diff
    Next X
End If
If F1MDF = 0 Then
    Form1.VScroll5.Top = Form1.Picture17.Top + Form1.Picture17.Height + 50
    Form1.VScroll5.Height = Form1.Picture29.Top + Form1.Picture29.Height - Form1.VScroll5.Top

    If Form1.VScroll5.Top > Form1.Picture18.Top + Form1.Picture18.Height + 50 Then
        Form1.SSPanel15.Top = Form1.VScroll5.Top
    Else
        Form1.SSPanel15.Top = Form1.Picture18.Top + Form1.Picture18.Height + 50
    End If

    Form1.SSPanel15.Height = Form1.Picture29.Top + Form1.Picture29.Height - Form1.SSPanel15.Top
End If
Dim LOSpace As Long
LOSpace = Form1.SSPanel15.Height - ((Form1.Label7(0).Height - 50) * 3)
LOSpace = LOSpace / 4
XX = Form1.Label7(0).Height
If LOSpace < 0 Then
    Extra = (LOSpace * 4) / 2
    LOSpace = 0
    
Else
    Extra = 0
End If
For X = 0 To 2
    'Form1.Label7(X).Top = 100 + ((Form1.SSPanel15.Height - 100) / 3) * (X)
    Form1.Label7(X).Top = ((Form1.Label7(0).Height - 50) * (X)) + (Extra * X) + LOSpace * (X + 1)
    
Next X

Call DoText(Form1.Picture18, Form1.Picture18.Font, Cap1, 270)

Form1.Picture17.ScaleMode = 3
Form1.Picture17.AutoRedraw = True

PH = Form1.Picture17.ScaleHeight
If MaxN >= 3 Or TF = 1 Or TF = 2 Or DiscF = 1 Then
    For X = 0 To PH
        Form1.Picture17.Line (0, X)-(Form1.Picture17.Width, X), HeatMap(CScale, 1020 - CInt((X / PH) * 1020))
        Form1.Picture26.Line (Form1.Picture26.ScaleWidth - 15, X)-(Form1.Picture26.ScaleWidth, X), HeatMap(CScale, 1020 - CInt((X / PH) * 1020))
    Next X
ElseIf MaxN = 2 Then
    For X = 0 To PH / 3
        Form1.Picture17.Line (0, X)-(Form1.Picture17.Width, X), HeatMap(CScale, 1020)
        Form1.Picture26.Line (Form1.Picture26.ScaleWidth - 15, X)-(Form1.Picture26.ScaleWidth, X), HeatMap(CScale, 1020)
    Next X
    For X = PH / 3 To (PH / 3) * 2
        Form1.Picture17.Line (0, X)-(Form1.Picture17.Width, X), HeatMap(CScale, 510)
        Form1.Picture26.Line (Form1.Picture26.ScaleWidth - 15, X)-(Form1.Picture26.ScaleWidth, X), HeatMap(CScale, 510)
    Next X
    For X = (PH / 3) * 2 To PH
        Form1.Picture17.Line (0, X)-(Form1.Picture17.Width, X), HeatMap(CScale, 0)
        Form1.Picture26.Line (Form1.Picture26.ScaleWidth - 15, X)-(Form1.Picture26.ScaleWidth, X), HeatMap(CScale, 0)
    Next X
ElseIf MaxN = 1 Then
    Form1.Picture17.Line (0, 0)-(Form1.Picture17.Width, PH / 2), HeatMap(CScale, 1020), BF
    Form1.Picture17.Line (0, PH / 2)-(Form1.Picture17.Width, PH), HeatMap(CScale, 0), BF
    Form1.Picture26.Line (Form1.Picture26.ScaleWidth - 15, 0)-(Form1.Picture26.ScaleWidth, PH / 2), HeatMap(CScale, 1020), BF
    Form1.Picture26.Line (Form1.Picture26.ScaleWidth - 15, PH / 2)-(Form1.Picture26.ScaleWidth, PH), HeatMap(CScale, 0), BF
End If
Form1.Picture17.ScaleMode = 1
Form1.Picture17.Refresh

For X = 0 To 4
    Form1.Label6(X).Visible = True
Next X
MaxN = tMaxN
End Sub

Public Function GetCritChi(Target, DF As Integer)
hPrb = Target * 1.000000001
lPrb = Target * 0.999999999
Dim HMChi As Double
'Calculate critical Chi
HMChi = 5
LastChi = 1
TPVal = 10
Do While TPVal > lPrb
   HMChi = HMChi * 2
   oHMChi = HMChi
   TPVal = chi2(HMChi, DF)
   HMChi = oHMChi
   If TPVal = 10 ^ -20 Then
      lPrb = TPVal
      Exit Do
    End If
    
Loop
LastLoChi = 0
LastHiChi = HMChi * 2
X = X
Do
    oHMChi = HMChi
    TPVal = chi2(HMChi, DF)
     HMChi = oHMChi
    If TPVal < lPrb Then
        TempChi = HMChi
        HMChi = HMChi - (HMChi - LastLoChi) / 2
        
        If HMChi = TempChi Then Exit Do
        LastHiChi = TempChi
    ElseIf TPVal > hPrb Then
        TempChi = HMChi
        HMChi = HMChi + (LastHiChi - HMChi) / 2
        oHMChi = HMChi
        If HMChi = TempChi Then Exit Do
        LastLoChi = TempChi
    Else
        Exit Do
    End If
    
Loop

'XX = chi2(1.5, 4)
X = X
GetCritChi = HMChi

End Function

Public Sub DrawRegionMatrix()
Dim RSize As Long, SConvert As Double
Dim UseAll As Byte
UseAll = 1
RSize = Len(StrainSeq(0))
If RSize > 1000 Then RSize = 1000
SConvert = RSize / Len(StrainSeq(0))





If DoneMatX(3) = 0 Then
    ReDim MatrixRR(RSize + 1, RSize + 1)
    Dim Size As Long, d As Long, P1 As Long, P2 As Long, Win As Long
        
    
    Dim PermutationX As Long
    PermutationX = 1000
    ENumb = SEventNumber
    
    Dim Excl() As Byte, Enu As Long, NC As Long
    ReDim Excl(SEventNumber), pCount(SEventNumber, AddNum)
        
    Dim BPV() As Double
    ReDim BPV(SEventNumber, AddNum)
    For X = 1 To SEventNumber
        For Y = 0 To AddNum
            BPV(X, Y) = LowestProb
        Next Y
    Next X
        
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            If XOverList(X, Y).ProgramFlag <= AddNum - 1 Then
                If BPV(SuperEventList(XOverList(X, Y).Eventnumber), XOverList(X, Y).ProgramFlag) > XOverList(X, Y).Probability And XOverList(X, Y).Probability > 0 Then
                    BPV(SuperEventList(XOverList(X, Y).Eventnumber), XOverList(X, Y).ProgramFlag) = XOverList(X, Y).Probability
                End If
            End If
        Next Y
    Next X
    
    Dim PValCon As Double
    For X = 1 To SEventNumber
        For Y = 0 To AddNum - 1
            If BPV(X, Y) = LowestProb And Confirm(X, Y) > 0 Then
                
                PValCon = ConfirmP(X, Y) / Confirm(X, Y)
                PValCon = 10 ^ (-PValCon)
                
                'ConfirmP(X, Y) = PValCon
                If BPV(X, Y) > PValCon Then
                    BPV(X, Y) = PValCon
                End If
            End If
        Next Y
    Next X
        
        
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            Enu = SuperEventList(XOverList(X, Y).Eventnumber)
                
            NC = 0
            If BPV(Enu, 0) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 1) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 3) < LowestProb Then
                    
                NC = NC + 1
                    
            ElseIf BPV(Enu, 4) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 2) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
            If BPV(Enu, 5) < LowestProb Then
                    
                NC = NC + 1
                    
            End If
                
            If XOverList(X, Y).Probability > 0 And XOverList(X, Y).Probability < LowestProb And NC > ConsensusProg And XOverList(X, Y).Accept <> 2 And XOverList(X, Y).ProgramFlag <= AddNum - 1 And XOverList(X, Y).MissIdentifyFlag <> 3 And XOverList(X, Y).MissIdentifyFlag <> 13 Then
                Excl(Enu) = 1
            
            End If
        Next Y
    Next X
        
        
    Dim BPos() As Long, APos() As Long, TypeSeq
    ReDim APos(Len(StrainSeq(0))), BPos(Len(StrainSeq(0)))
    
    TypeSeq = 0
    If TypeSeq > PermNextNo Then TypeSeq = 0
        
    For X = 1 To Len(StrainSeq(0))
        APos(X) = X - SeqSpaces(X, TypeSeq)
        BPos(X - SeqSpaces(X, TypeSeq)) = X
    Next X
        
    If X = 12345 Then
        Win = 200
        B = 0
        C = 0
        If BPCvalFlag = 0 Or X = X Then
            BPCvalFlag = 1
            ReDim BPCVal(1, 1)
            Dim DN As Long
            DN = 1
            Call RecombMapPermsB(DN, APos(), BPos(), Excl(), BPCVal(), Win, PermutationX, 3)
        End If
    End If
    
    zz = 0
        
    'Make exclusions
    Dim RS1 As Long, LS1 As Long, RS2 As Long, LS2 As Long
    SS = GetTickCount
    Dim LastTick As Long, CurTick As Long
    LastTick = 0
    For X = 1 To ENumb
        
        If Excl(X) = 1 And (BestEvent(X, 0) > 0 Or BestEvent(X, 1) > 0) Then
            zz = zz + 1
            CNum = 0
            sen = X 'SuperEventlist(XOverList(BestEvent(X, 0), BestEvent(X, 1)).Eventnumber)
               
            If XOverList(BestEvent(X, 0), BestEvent(X, 1)).Accept = 1 Or (UseAll = 1 And XOverList(BestEvent(X, 0), BestEvent(X, 1)).Accept <> 2) Then
                d = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Daughter
                P1 = XOverList(BestEvent(X, 0), BestEvent(X, 1)).MajorP
                P2 = XOverList(BestEvent(X, 0), BestEvent(X, 1)).MinorP
                ST = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Beginning
                EN = XOverList(BestEvent(X, 0), BestEvent(X, 1)).Ending
                
                If X = X Then
                    Dummy = MakeRecCMatrix(Len(StrainSeq(0)), RSize, ST, EN, SConvert, MatrixRR(0, 0))
                Else
                    LS1 = -1: LS2 = -1
                    If ST < EN Then
                        For A = 1 To ST - 1
                            RS1 = CInt(A * SConvert)
                            If LS1 <> RS1 Then
                                LS1 = RS1
                            
                                For B = ST To EN
                                    RS2 = CInt(B * SConvert)
                                    If LS2 <> RS2 Then
                                        LS2 = RS2
                                        MatrixRR(RS2, RS1) = MatrixRR(RS2, RS1) + 1
                                        MatrixRR(RS1, RS2) = MatrixRR(RS1, RS2) + 1
                                    End If
                                Next B
                            End If
                        Next A
                        For A = EN + 1 To Len(StrainSeq(0))
                            RS1 = CInt(A * SConvert)
                            If LS1 <> RS1 Then
                                LS1 = RS1
                            
                                For B = ST To EN
                                    RS2 = CInt(B * SConvert)
                                    If LS2 <> RS2 Then
                                        LS2 = RS2
                                        MatrixRR(RS2, RS1) = MatrixRR(RS2, RS1) + 1
                                        MatrixRR(RS1, RS2) = MatrixRR(RS1, RS2) + 1
                                    End If
                                Next B
                            End If
                        Next A
                    Else
                        For A = EN + 1 To ST - 1
                            RS1 = CInt(A * SConvert)
                            If LS1 <> RS1 Then
                                LS1 = RS1
                            
                                For B = 1 To EN
                                    RS2 = CInt(B * SConvert)
                                    If LS2 <> RS2 Then
                                        LS2 = RS2
                                        MatrixRR(RS2, RS1) = MatrixRR(RS2, RS1) + 1
                                        MatrixRR(RS1, RS2) = MatrixRR(RS1, RS2) + 1
                                    End If
                                Next B
                                For B = ST To Len(StrainSeq(0))
                                    RS2 = CInt(B * SConvert)
                                    If LS2 <> RS2 Then
                                        LS2 = RS2
                                        MatrixRR(RS2, RS1) = MatrixRR(RS2, RS1) + 1
                                        MatrixRR(RS1, RS2) = MatrixRR(RS1, RS2) + 1
                                    End If
                                Next B
                            End If
                        Next A
                        
                    
                    End If
                End If
            End If
        End If
        CurTick = GetTickCount
        If Abs(CurTick - LastTick) > 200 Then
            LastTick = CurTick
            Form1.SSPanel1.Caption = Trim(Str(X)) + " of " + Trim(Str(SEventNumber)) + " events mapped"
            Form1.ProgressBar1 = (X / SEventNumber) * 100
        End If
    Next X
    Form1.ProgressBar1 = 100
    EE = GetTickCount
    TT = EE - SS '2.687,2.704
    '3.77,2.766, 2.797
    MaxN = FindMaxN(RSize, MatrixRR(0, 0))
    MatBound(3) = MaxN
    Form1.SSPanel1.Caption = "Drawing matrix"
Else
    MaxN = MatBound(3)
End If
SS = GetTickCount


    Form1.Picture26.ScaleMode = 3
    Form1.Picture26.AutoRedraw = True
    Dim XAddj As Double, PosS(1) As Long, PosE(1) As Long, DistD As Long
        
        DistD = RSize / MatZoom(3)
        PosS(0) = MatCoord(3, 0)
        PosE(0) = PosS(0) + DistD
        PosS(1) = MatCoord(3, 1)
        PosE(1) = PosS(1) + DistD
        If PosE(1) > (UBound(MatrixRR, 1) - 1) - 1 Then PosE(1) = (UBound(MatrixRR, 1) - 1) - 1
        If PosE(0) > (UBound(MatrixRR, 1) - 1) - 1 Then PosE(0) = (UBound(MatrixRR, 1) - 1) - 1
       
        If PosS(0) < 0 Then SX = 0 Else SX = PosS(0)
        If PosS(1) < 0 Then SY = 0 Else SY = PosS(1)
        
        If MaxN > 2 Then
            MaxN = CLng((MaxN / 4) + 0.5) * 4
        Else
            MaxN = 2
        End If
        Form1.Picture26.Picture = LoadPicture()
        XAddj = (Form1.Picture26.ScaleHeight) / DistD
        Pict = Form1.Picture26.hDC
        
        If X = 12345 Then
            MakeHeatPlot MaxN, UBound(HeatMap, 1), CurScale, Pict, RSize, XAddj, HeatMap(0, 0), MatrixRR(0, 0)
        Else
            If X = X Then
                Call DrawmatsVB(0, PosE(), PosS(), SX, SY, XAddj, MatrixRR(), HeatMap(), CurScale, MaxN)
            Else
                RangeN = MaxN - MinN
                For X = 1 To RSize 'Step 20
                    For Y = X To RSize 'Step 20
                        SetPixelV Pict, CInt(X * XAddj), CInt(Y * XAddj), HeatMap(CurScale, CInt((MatrixRR(X, Y) / MaxN) * 1020))
                        SetPixelV Pict, CInt(Y * XAddj), CInt(X * XAddj), HeatMap(CurScale, CInt((MatrixRR(X, Y) / MaxN) * 1020))
                    Next Y
                    Form1.Picture26.Refresh
                    Form1.ProgressBar1 = 70 + (X / RSize) * 30
                Next X
            End If
        End If
        If DontDoKey = 0 Then
            Call DoKey(0, MaxN, MinN, 3, "Number of events", CurScale)
        End If
        DoneMatX(3) = 1
        
        Form1.Picture26.Refresh
    Form1.ProgressBar1 = 0
    Form1.SSPanel1 = ""
End Sub


Public Sub CheckDrop(Steps() As Long, SEventNumber As Long, StepNo As Long, Nextno As Long, oNextno As Long, NumRecsI() As Long, RedoListSize As Long, RedoList() As Long)

Dim X As Long, z As Long, A As Long

X = Nextno
Do While X > oNextno
            
    If NumRecsI(X) = 0 Then
        If RedoListSize > 0 Then
            z = 0
            Do While z <= RedoListSize
                For A = 1 To 3
                    If RedoList(A, z) = X Then
                        Exit For
                    End If
                Next A
                If A < 4 Then
                    For A = 0 To 3
                        RedoList(A, z) = RedoList(A, RedoListSize)
                    Next A
                    RedoListSize = RedoListSize - 1
                Else
                    z = z + 1
                End If
            Loop
        End If
                    
        Steps(0, StepNo) = 3 'ie delete a sequence ....
        Steps(1, StepNo) = Nextno 'this is the seqence.....
        Steps(4, StepNo) = SEventNumber + 1
        StepNo = StepNo + 1
        UB = UBound(Steps, 2)
        If StepNo > UB Then
            ReDim Preserve Steps(4, UB + 100)
        End If
                    
        Nextno = Nextno - 1
        X = X - 1
    Else
        Exit Do
    End If
                
Loop
End Sub

Public Sub CalcMatchX(ISeqs() As Long, CompMat() As Long, OKSeq() As Double, BPos As Long, EPos As Long)
Dim VarSiteSmooth() As Double, LenVarSeq As Long, VarSiteMap() As Byte, VXPos() As Long, SWin As Long, VarSiteSmooth2() As Double, Tot As Double
ReDim VarSiteMap(2, Len(StrainSeq(0)), Nextno), VXPos(Len(StrainSeq(0)))
LenVarSeq = 0
Dim VRPos() As Long
ReDim VRPos(Len(StrainSeq(0)) + 1)
VRPos(Len(StrainSeq(0)) + 1) = Len(StrainSeq(0)) + 1

Dim CntHit() As Double, NCnt As Long
ReDim CntHit(2, 1, Nextno)
ReDim VarSiteSmooth(2, Len(StrainSeq(0)), Nextno)
'ReDim VarSiteSmoothxx(2, Len(StrainSeq(0)), NextNo)
If X = X Then
    LenVarSeq = MakeVarMap(Nextno, Len(StrainSeq(0)), SeqNum(0, 0), VarSiteMap(0, 0, 0), VRPos(0), VXPos(0), ISeqs(0), CompMat(0, 0))
Else
    For X = 0 To Len(StrainSeq(0))
        VRPos(X) = LenVarSeq
        If SeqNum(X, ISeqs(0)) <> 46 Then
            If SeqNum(X, ISeqs(1)) <> 46 Then
                If SeqNum(X, ISeqs(2)) <> 46 Then
                    If SeqNum(X, ISeqs(0)) <> SeqNum(X, ISeqs(1)) Or SeqNum(X, ISeqs(0)) <> SeqNum(X, ISeqs(2)) Then
                        LenVarSeq = LenVarSeq + 1
                        VXPos(LenVarSeq) = X
                        For Y = 0 To 2
                            For z = 0 To Nextno
                                If SeqNum(X, ISeqs(Y)) = SeqNum(X, z) Then
                                    VarSiteMap(Y, LenVarSeq, z) = 2
                                ElseIf SeqNum(X, z) <> SeqNum(X, ISeqs(CompMat(Y, 0))) And SeqNum(X, z) <> SeqNum(X, ISeqs(CompMat(Y, 1))) Then
                                    VarSiteMap(Y, LenVarSeq, z) = 1
                                End If
                            Next z        '    VarSiteMap(LenVarSeq) = 1
                        Next Y                'End If
                    End If
                End If
            End If
        End If
    Next X
End If



'smooth varsitemap
SWin = 10
If X = X Then
    Dummy = MakeCntHit(BPos, EPos, SWin, Nextno, LenVarSeq, Len(StrainSeq(0)), CntHit(0, 0, 0), VarSiteMap(0, 0, 0), VarSiteSmooth(0, 0, 0), VRPos(0))
    
Else
    For A = 0 To 2
        For B = 0 To Nextno
            Tot = 0
            For X = 1 - SWin To 1 + SWin
            
                If X < 1 Then
                    z = (LenVarSeq + X) '2450
                ElseIf X > LenVarSeq Then
                    z = X - LenVarSeq
                Else
                    z = X
                End If
                Tot = Tot + VarSiteMap(A, z, B)
            Next X
            
            VarSiteSmooth(A, 1, B) = Tot / ((SWin * 2 + 1) * 2)
            For X = 2 To LenVarSeq
                z = X - SWin - 1
                If z < 1 Then
                    Tot = Tot - VarSiteMap(A, (LenVarSeq + z), B)
                ElseIf z > LenVarSeq Then
                    Tot = Tot - VarSiteMap(A, z - LenVarSeq, B)
                Else
                    Tot = Tot - VarSiteMap(A, z, B)
                End If
                z = X + SWin
                If z > LenVarSeq Then
                    Tot = Tot + VarSiteMap(A, z - LenVarSeq, B)
                Else
                    Tot = Tot + VarSiteMap(A, z, B)
                End If
                VarSiteSmooth(A, X, B) = Tot / ((SWin * 2 + 1) * 2)
                
            Next X
        Next B
    Next A

'End If

    If BPos > 1 Then
        ST = BPos - 1
    Else
        ST = BPos
    End If
    
    
    NCnt = 0
    If BPos < EPos Then
        For X = 1 To VRPos(ST - 1)
            NCnt = NCnt + 1
            For Y = 0 To 2
                For z = 0 To Nextno
                    If VarSiteSmooth(Y, X, z) > 0.6 Then
                       CntHit(Y, 0, z) = CntHit(Y, 0, z) + (VarSiteSmooth(Y, X, z) - 0.6) / 0.4
                       X = X
                    End If
                Next z
            Next Y
            
        Next X
        For X = VRPos(EPos + 1) To VRPos(Len(StrainSeq(0)))
            NCnt = NCnt + 1
            For Y = 0 To 2
                For z = 0 To Nextno
                    If VarSiteSmooth(Y, X, z) > 0.6 Then
                        CntHit(Y, 0, z) = CntHit(Y, 0, z) + (VarSiteSmooth(Y, X, z) - 0.6) / 0.4
                    End If
                Next z
            Next Y
        Next X
        
        For Y = 0 To 2
            For z = 0 To Nextno
            'If Z = 30 And Y = 1 Then
            '       X = X
            '    End If
                CntHit(Y, 0, z) = CntHit(Y, 0, z) / NCnt
                X = X
            Next z
        Next Y
        NCnt = 0
        For X = VRPos(BPos) To VRPos(EPos)
            NCnt = NCnt + 1
            For Y = 0 To 2
                For z = 0 To Nextno
                    If VarSiteSmooth(Y, X, z) > 0.6 Then
                        CntHit(Y, 1, z) = CntHit(Y, 1, z) + (VarSiteSmooth(Y, X, z) - 0.6) / 0.4
                    End If
                Next z
            Next Y
        Next X
        For Y = 0 To 2
            For z = 0 To Nextno
              '  If Z = 30 And Y = 1 Then
              '     X = X
              '  End If
                CntHit(Y, 1, z) = CntHit(Y, 1, z) / NCnt
                X = X
            Next z
        Next Y
    Else
        For X = 1 To VRPos(EPos)
            NCnt = NCnt + 1
            For Y = 0 To 2
                For z = 0 To Nextno
                    If VarSiteSmooth(Y, X, z) > 0.6 Then
                       CntHit(Y, 1, z) = CntHit(Y, 1, z) + (VarSiteSmooth(Y, X, z) - 0.6) / 0.4
                    End If
                Next z
            Next Y
            
        Next X
        For X = VRPos(BPos) To VRPos(Len(StrainSeq(0)))
            NCnt = NCnt + 1
            For Y = 0 To 2
                For z = 0 To Nextno
                    If VarSiteSmooth(Y, X, z) > 0.6 Then
                        CntHit(Y, 1, z) = CntHit(Y, 1, z) + (VarSiteSmooth(Y, X, z) - 0.6) / 0.4
                    End If
                Next z
            Next Y
        Next X
        For Y = 0 To 2
            For z = 0 To Nextno
             '   If Z = 30 And Y = 1 Then
             '      X = X
             '   End If
                CntHit(Y, 1, z) = CntHit(Y, 1, z) / NCnt
                X = X
            Next z
        Next Y
        NCnt = 0
        For X = VRPos(EPos + 1) To VRPos(BPos - 1)
            NCnt = NCnt + 1
            For Y = 0 To 2
                For z = 0 To Nextno
                    If VarSiteSmooth(Y, X, z) > 0.6 Then
                        CntHit(Y, 0, z) = CntHit(Y, 0, z) + (VarSiteSmooth(Y, X, z) - 0.6) / 0.4
                    End If
                Next z
            Next Y
        Next X
        For Y = 0 To 2
            For z = 0 To Nextno
                'If Z = 30 And Y = 1 Then
                '   X = X
                'End If
                CntHit(Y, 0, z) = CntHit(Y, 0, z) / NCnt
                    
            Next z
        Next Y
    End If
    
    
End If

For Y = 0 To 2
    For z = 0 To Nextno
       
        OKSeq(Y, 17, z) = CntHit(Y, 0, z) * CntHit(Y, 1, z)
           
    Next z
Next Y

'Call SmoothIt(VarsiteSmooth(), VarSiteSmooth2(), 3, LenVarSeq, 0.25)
'Form1.Picture7.AutoRedraw = True

End Sub

Public Sub DrawHorstInterval(ST, LP)

HorstFlag = 1
Dim Score As Double, ScoreTable() As Double
ReDim ScoreTable(255, 255, 255)
ScoreTable(66, 66, 66) = 1.92
ScoreTable(66, 66, 68) = 2.42
ScoreTable(66, 66, 72) = 2.42
ScoreTable(66, 66, 85) = 1.25


ScoreTable(66, 68, 66) = 3.46
ScoreTable(66, 68, 68) = 3.94
ScoreTable(66, 68, 72) = 4.3
ScoreTable(66, 68, 85) = 3.46


ScoreTable(66, 72, 66) = 3.46
ScoreTable(66, 72, 68) = 4.3
ScoreTable(66, 72, 72) = 3.94
ScoreTable(66, 72, 85) = 3.46

ScoreTable(66, 85, 66) = 0.04
ScoreTable(66, 85, 68) = 1.13
ScoreTable(66, 85, 72) = 1.13
ScoreTable(66, 85, 85) = 0.61


ScoreTable(68, 66, 66) = 2.95
ScoreTable(68, 66, 68) = 3.46
ScoreTable(68, 66, 72) = 3.46
ScoreTable(68, 66, 85) = 2.3

ScoreTable(68, 68, 66) = 4.43
ScoreTable(68, 68, 68) = 4.91
ScoreTable(68, 68, 72) = 5.27
ScoreTable(68, 68, 85) = 4.43

ScoreTable(68, 72, 66) = 5.14
ScoreTable(68, 72, 68) = 5.98
ScoreTable(68, 72, 72) = 5.63
ScoreTable(68, 72, 85) = 5.14

ScoreTable(68, 85, 66) = 2.3
ScoreTable(68, 85, 68) = 3.46
ScoreTable(68, 85, 72) = 3.46
ScoreTable(68, 85, 85) = 2.95


ScoreTable(72, 66, 66) = 2.95
ScoreTable(72, 66, 68) = 3.46
ScoreTable(72, 66, 72) = 3.46
ScoreTable(72, 66, 85) = 2.3

ScoreTable(72, 68, 66) = 5.14
ScoreTable(72, 68, 68) = 5.63
ScoreTable(72, 68, 72) = 5.98
ScoreTable(72, 68, 85) = 5.14

ScoreTable(72, 72, 66) = 4.43
ScoreTable(72, 72, 68) = 5.27
ScoreTable(72, 72, 72) = 4.91
ScoreTable(72, 72, 85) = 4.43

ScoreTable(72, 85, 66) = 2.3
ScoreTable(72, 85, 68) = 3.46
ScoreTable(72, 85, 72) = 3.46
ScoreTable(72, 85, 85) = 2.95

ScoreTable(85, 66, 66) = 0.61
ScoreTable(85, 66, 68) = 1.13
ScoreTable(85, 66, 72) = 1.13
ScoreTable(85, 66, 85) = 0.04

ScoreTable(85, 68, 66) = 3.46
ScoreTable(85, 68, 68) = 3.94
ScoreTable(85, 68, 72) = 4.3
ScoreTable(85, 68, 85) = 3.46

ScoreTable(85, 72, 66) = 3.46
ScoreTable(85, 72, 68) = 4.3
ScoreTable(85, 72, 72) = 3.94
ScoreTable(85, 72, 85) = 3.46

ScoreTable(85, 85, 66) = 1.25
ScoreTable(85, 85, 68) = 2.42
ScoreTable(85, 85, 72) = 2.42
ScoreTable(85, 85, 85) = 1.92


Score = 0
Dim ScorePlot() As Double, Pos(2)
ReDim ScorePlot(2, Len(StrainSeq(0)))
For z = 0 To 2
    For X = ST + z To ST + LP - 6 Step 3
        'XX = SeqNum(X, 0) ' ='a=66,c=68, g=72, t=85
        Y = 0
        Pos(0) = 0: Pos(1) = 0: Pos(2) = 0
        Do
            If SeqNum(z + Y, 0) <> 46 Then
                Pos(Y) = X + Y
                Y = Y + 1
                If Y = 3 Then Exit Do
                If Y + X > ST + LP Then Exit For
                If Y + X > Len(StrainSeq(0)) Then Exit For
            End If
        Loop
        ScorePlot(z, Pos(0)) = ScoreTable(SeqNum(Pos(0), 0), SeqNum(Pos(1), 0), SeqNum(Pos(2), 0))
        ScorePlot(z, Pos(1)) = ScoreTable(SeqNum(Pos(0), 0), SeqNum(Pos(1), 0), SeqNum(Pos(2), 0))
        ScorePlot(z, Pos(2)) = ScoreTable(SeqNum(Pos(0), 0), SeqNum(Pos(1), 0), SeqNum(Pos(2), 0))
    Next X
Next z



    Dim YScaleFactor As Double

    Form1.Picture7.Picture = LoadPicture()
    
    YScaleFactor = 0.85
    PicHeight = Form1.Picture7.Height * YScaleFactor
    Form1.Picture7.Cls
    'Draw homology plot in picturebox 7

    Dim PntAPI As POINTAPI
    Dim Pict As Long

    
    XFactor = ((Form1.Picture7.Width - 40) / (LP))
    Call DoAxes(0, ST + LP, -1, 6, 0, 0, "Energy")
    XFactor = ((Form1.Picture7.Width - 40) / (LP))
    Pict = Form1.Picture7.hDC
    For z = 0 To 2
        If z = 0 Then
            Form1.Picture7.ForeColor = RGB(255, 0, 0)
        ElseIf z = 1 Then
            Form1.Picture7.ForeColor = RGB(0, 255, 0)
        ElseIf z = 2 Then
            Form1.Picture7.ForeColor = RGB(0, 0, 255)
        End If
        MoveToEx Pict, 30 + 1 * XFactor, PicHeight - (15 + (ScorePlot(z, 1) / 6) * (PicHeight - 35)), PntAPI

        For X = ST + 1 To ST + LP - 1 'Len(StrainSeq(0))
            LineTo Pict, 30 + (X - (ST - 1)) * XFactor, PicHeight - (15 + (ScorePlot(z, X) / 6) * (PicHeight - 35))
        Next X
    Next z
    Form1.Picture7.Refresh
End Sub


Public Sub CalcMatch(SX, SY, sa, SB)
Dim LenVarSeq As Long, VarSiteMap() As Byte, VXPos() As Long, SWin As Long, VarSiteSmooth2() As Double, VarSiteSmooth() As Double, Tot As Double
ReDim VarSiteMap(Len(StrainSeq(0))), VXPos(Len(StrainSeq(0)))
LenVarSeq = 0
oSA = sa
oSB = SB
OSX = SX
oSY = SY
sa = TreeTrace(sa)
SX = TreeTrace(SX)
SY = TreeTrace(SY)
SB = TreeTrace(SB)
'Call ModSeqNum(1)
For X = 0 To Len(StrainSeq(0))
            If SeqNum(X, sa) <> 46 Then
                If SeqNum(X, SB) <> 46 Then
                    If SeqNum(X, SX) <> 46 Then
                        'If SeqNum(X, SY) <> 46 Then
                            If SeqNum(X, sa) <> SeqNum(X, SB) Or SeqNum(X, sa) <> SeqNum(X, SX) Then
                            'If SeqNum(X, SA) <> SeqNum(X, SB) Or SeqNum(X, SA) <> SeqNum(X, SX) Or SeqNum(X, SA) <> SeqNum(X, SY) Then
                                'LenVarSeq = LenVarSeq + 1
                                'VXPos(LenVarSeq) = X
                                LenVarSeq = LenVarSeq + 1
                                VXPos(LenVarSeq) = X
                                If SeqNum(X, SX) = SeqNum(X, SY) Then
                                    VarSiteMap(LenVarSeq) = 2
                                ElseIf SeqNum(X, SY) <> SeqNum(X, SB) And SeqNum(X, SY) <> SeqNum(X, sa) Then
                                    VarSiteMap(LenVarSeq) = 1
                                End If
                                'If SeqNum(X, SA) <> SeqNum(X, SB) Then
                                '    LenVarSeq = LenVarSeq + 1
                                '    VXPos(LenVarSeq) = X
                                '    If SeqNum(X, SX) = SeqNum(X, SY) Then
                                '        VarSiteMap(LenVarSeq) = 1
                                '    End If
                                    'If SeqNum(X, SX) <> SeqNum(X, SY) Then
                                    '
                                    '    If SeqNum(X, SX) = SeqNum(X, SA) And SeqNum(X, SY) = SeqNum(X, SB) Then
                                    '        VarSiteMap(LenVarSeq) = 1
                                    '
                                    '    ElseIf SeqNum(X, SX) = SeqNum(X, SB) And SeqNum(X, SY) = SeqNum(X, SA) Then
                                    '        VarSiteMap(LenVarSeq) = 1
                                    '        'LenVarSeq = LenVarSeq + 1
                                    '        'VXPos(LenVarSeq) = X
                                    '    'ElseIf SeqNum(X, SX) = SeqNum(X, SA) Or SeqNum(X, SX) = SeqNum(X, SB) Or SeqNum(X, SY) = SeqNum(X, SA) Or SeqNum(X, SY) = SeqNum(X, SB) Then
                                    '    '    VarSiteMap(LenVarSeq) = 1
                                    '
                                    '    End If
                                    'End If
                                'If SeqNum(X, SX) <> SeqNum(X, SY) Then
                                '    VarSiteMap(LenVarSeq) = 1
                                'End If
                            End If
                        'End If
                    End If
                End If
            End If
Next X
'Call UnModSeqNum(1)
ReDim VarSiteSmooth(Len(StrainSeq(0)))
ReDim VarSiteSmooth2(Len(StrainSeq(0)))
'smooth varsitemap
SWin = 10
Tot = 0
XX = 0
For X = 1 - SWin To 1 + SWin

    If X < 1 Then
        z = (LenVarSeq + X) '2450
    ElseIf X > LenVarSeq Then
        z = X - LenVarSeq
    Else
        z = X
    End If
    Tot = Tot + VarSiteMap(z)
Next X

VarSiteSmooth(1) = Tot / ((SWin * 2 + 1) * 2)
For X = 2 To LenVarSeq
    z = X - SWin - 1
    If z < 1 Then
        Tot = Tot - VarSiteMap((LenVarSeq + z))
    ElseIf z > LenVarSeq Then
        Tot = Tot - VarSiteMap(z - LenVarSeq)
    Else
        Tot = Tot - VarSiteMap(z)
    End If
    z = X + SWin
    If z > LenVarSeq Then
        Tot = Tot + VarSiteMap(z - LenVarSeq)
    Else
        Tot = Tot + VarSiteMap(z)
    End If
    VarSiteSmooth(X) = Tot / ((SWin * 2 + 1) * 2)
    X = X
Next X

Call SmoothIt(VarSiteSmooth(), VarSiteSmooth2(), 3, LenVarSeq, 0.25)
'Form1.Picture7.AutoRedraw = True
YScaleFactor = 0.85
    PicHeight = Form1.Picture7.Height * YScaleFactor
    
    XFactor = ((Form1.Picture7.Width - 40) / Len(StrainSeq(0)))

sa = oSA
SB = oSB
SX = OSX
SY = oSY

Dim Pict As Long, PntAPI As POINTAPI
' MoveToEx pict, 30 + XDiffPos(1) * XFactor, PicHeight - (15 + (((XOverHomologyNum(LenXOverSeq, Y) - RDPLD) / (RDPUD - RDPLD))) * (PicHeight - 35)), PntAPI
Form1.Picture7.ForeColor = RGB(170, 170, 170)
Form1.Picture7.DrawWidth = 2
Pict = Form1.Picture7.hDC
MoveToEx Pict, 30 + VXPos(1) * XFactor, 8, PntAPI
LineTo Pict, 30 + VXPos(1) * XFactor, 9
'MoveToEx pict, 30 + VXPos(1) * XFactor, PicHeight - (15 + VarSiteSmooth2(1) * (PicHeight - 35)), PntAPI
For X = 2 To LenVarSeq
        If X = X Then
            If (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.55 Then
                Form1.Picture7.ForeColor = HeatMap(0, 1020 - CLng((((VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 - 0.55) / 0.45) * 1020))
            Else
                Form1.Picture7.ForeColor = HeatMap(0, 1020)
            End If
            Pict = Form1.Picture7.hDC
        Else
            If (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.95 Then
                Form1.Picture7.ForeColor = RGB(0, 0, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.9375 Then
                Form1.Picture7.ForeColor = RGB(0, 32, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.925 Then
                Form1.Picture7.ForeColor = RGB(0, 64, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.9125 Then
                Form1.Picture7.ForeColor = RGB(0, 92, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.9 Then
                Form1.Picture7.ForeColor = RGB(0, 128, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.875 Then
                Form1.Picture7.ForeColor = RGB(0, 160, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.875 Then
                Form1.Picture7.ForeColor = RGB(0, 192, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.8625 Then
                Form1.Picture7.ForeColor = RGB(0, 224, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.85 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 255)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.8375 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 224)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.825 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 198)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.8125 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 160)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.8 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 128)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.7875 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 96)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.775 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 64)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.7625 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 32)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.75 Then
                Form1.Picture7.ForeColor = RGB(0, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.7375 Then
                Form1.Picture7.ForeColor = RGB(32, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.725 Then
                Form1.Picture7.ForeColor = RGB(64, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.7125 Then
                Form1.Picture7.ForeColor = RGB(96, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.7 Then
                Form1.Picture7.ForeColor = RGB(128, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.6875 Then
                Form1.Picture7.ForeColor = RGB(160, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.675 Then
                Form1.Picture7.ForeColor = RGB(192, 255, 0)
                Pict = Form1.Picture7.hDC
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.6625 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(224, 255, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.65 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 255, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.6375 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 224, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.625 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 192, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.6125 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 160, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.6 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 128, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.5875 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 96, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.575 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 64, 0)
            ElseIf (VarSiteSmooth2(X) + VarSiteSmooth2(X - 1)) / 2 > 0.5625 Then
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 32, 0)
            Else
                Pict = Form1.Picture7.hDC
                Form1.Picture7.ForeColor = RGB(255, 0, 0)
            End If
        End If
        LineTo Pict, 30 + VXPos(X) * XFactor, 8
        LineTo Pict, 30 + VXPos(X) * XFactor, 9
        'LineTo pict, (30 + VXPos(X) * XFactor), (PicHeight - (15 + VarSiteSmooth2(X) * (PicHeight - 35)))
Next X
Form1.Picture7.DrawWidth = 1
Form1.Picture7.Refresh
End Sub
Public Sub SmoothIt(VarSiteSmooth() As Double, VarSiteSmooth2() As Double, SWin, LenVarSeq, MinNo)
Dim Tot As Double
Tot = 0
For X = 1 - SWin To 1 + SWin

    If X < 1 Then
        z = (LenVarSeq + X) '2450
    ElseIf X > LenVarSeq Then
        z = X - LenVarSeq
    Else
        z = X
    End If
    Tot = Tot + VarSiteSmooth(z)
Next X

VarSiteSmooth2(1) = ((Tot / ((SWin * 2 + 1))) - MinNo) / (1 - MinNo)
If VarSiteSmooth2(1) < 0 Then VarSiteSmooth2(1) = 0
For X = 2 To LenVarSeq
    z = X - SWin - 1
    If z < 1 Then
        Tot = Tot - VarSiteSmooth((LenVarSeq + z))
    ElseIf z > LenVarSeq Then
        Tot = Tot - VarSiteSmooth(z - LenVarSeq)
    Else
        Tot = Tot - VarSiteSmooth(z)
    End If
    z = X + SWin
    If z > LenVarSeq Then
        Tot = Tot + VarSiteSmooth(z - LenVarSeq)
    Else
        Tot = Tot + VarSiteSmooth(z)
    End If
    VarSiteSmooth2(X) = ((Tot / ((SWin * 2 + 1))) - MinNo) / (1 - MinNo)
    If VarSiteSmooth2(X) < 0 Then VarSiteSmooth2(X) = 0
    
Next X
End Sub
Public Sub SmoothIt2D(VarSiteSmooth() As Double, VarSiteSmooth2() As Double, SWin, LenVarSeq, MinNo, numE)
Dim Tot As Double


For Y = 0 To numE
    Tot = 0
    For X = 1 - SWin To 1 + SWin
    
        If X < 1 Then
            z = (LenVarSeq + X) '2450
        ElseIf X > LenVarSeq Then
            z = X - LenVarSeq
        Else
            z = X
        End If
        Tot = Tot + VarSiteSmooth(Y, z)
    Next X
    
    VarSiteSmooth2(Y, 1) = ((Tot / ((SWin * 2 + 1))) - MinNo) / (1 - MinNo)
    If VarSiteSmooth2(Y, 1) < 0 Then VarSiteSmooth2(Y, 1) = 0
    For X = 2 To LenVarSeq
        z = X - SWin - 1
        If z < 1 Then
            Tot = Tot - VarSiteSmooth(Y, (LenVarSeq + z))
        ElseIf z > LenVarSeq Then
            Tot = Tot - VarSiteSmooth(Y, z - LenVarSeq)
        Else
            Tot = Tot - VarSiteSmooth(Y, z)
        End If
        z = X + SWin
        If z > LenVarSeq Then
            Tot = Tot + VarSiteSmooth(Y, z - LenVarSeq)
        Else
            Tot = Tot + VarSiteSmooth(Y, z)
        End If
        VarSiteSmooth2(Y, X) = ((Tot / ((SWin * 2 + 1))) - MinNo) / (1 - MinNo)
        If VarSiteSmooth2(Y, X) < 0 Then VarSiteSmooth2(Y, X) = 0
        
    Next X
Next Y
End Sub
Public Sub RCheckWithOther(S1, S2, S3)
        oSeq1 = Seq1
        oSeq2 = Seq2
        oSeq3 = Seq3
        Seq1 = S1
        Seq2 = S2
        Seq3 = S3
'If DontDoComboFlag = 0 Then
        Dim tSeqNum() As Integer
        If LongWindedFlag = 1 Then
           Call ModSeqNum(0)
           
        End If
        
        
        Form1.Refresh
        Form1.Picture7.ScaleMode = 3

        If XOverList(RelX, RelY).ProgramFlag = 8 Or XOverList(RelX, RelY).ProgramFlag = 8 + AddNum Or XOverList(RelX, RelY).ProgramFlag = 6 Or XOverList(RelX, RelY).ProgramFlag = 6 + AddNum Or ((XOverList(RelX, RelY).ProgramFlag = 1 Or XOverList(RelX, RelY).ProgramFlag = 1 + AddNum) And pGCTripletflag = 1) Or XOverList(RelX, RelY).ProgramFlag = 4 + AddNum Or XOverList(RelX, RelY).ProgramFlag = 4 Or XOverList(RelX, RelY).ProgramFlag = 5 Or XOverList(RelX, RelY).ProgramFlag = 5 + AddNum Or XOverList(RelX, RelY).ProgramFlag = 0 Or XOverList(RelX, RelY).ProgramFlag = 0 + AddNum Or XOverList(RelX, RelY).ProgramFlag = 2 Or XOverList(RelX, RelY).ProgramFlag = 2 + AddNum Or ((XOverList(RelX, RelY).ProgramFlag = 3 Or XOverList(RelX, RelY).ProgramFlag = 3 + AddNum) And MCTripletFlag = 0) Then
            If Form1.Combo1.ListIndex < 5 Then
                CurrentCheck = Form1.Combo1.ListIndex
            ElseIf Form1.Combo1.ListIndex = 5 Then
                CurrentCheck = 10
            ElseIf Form1.Combo1.ListIndex = 6 Then
                CurrentCheck = 5
            ElseIf Form1.Combo1.ListIndex = 7 Then
                CurrentCheck = 13
            ElseIf Form1.Combo1.ListIndex = 11 Then
                CurrentCheck = 8
            ElseIf Form1.Combo1.ListIndex = 12 Then
                CurrentCheck = 11
            ElseIf Form1.Combo1.ListIndex = 13 Then
                CurrentCheck = 12
            ElseIf Form1.Combo1.ListIndex = 14 Then
                CurrentCheck = 14
            ElseIf Form1.Combo1.ListIndex = 15 Then
                CurrentCheck = 15
            ElseIf Form1.Combo1.ListIndex = 9 Then
                CurrentCheck = 16
            Else
                CurrentCheck = Form1.Combo1.ListIndex - 2
            End If
        Else 'If CurrentCheck > -1 Then
            If Form1.Combo1.ListIndex = 0 Then
                CurrentCheck = 20
            ElseIf Form1.Combo1.ListIndex = 1 Then
                CurrentCheck = 21
            ElseIf Form1.Combo1.ListIndex = 2 Then
                CurrentCheck = 25
            ElseIf Form1.Combo1.ListIndex = 3 Then
                CurrentCheck = 22
                
            ElseIf Form1.Combo1.ListIndex = 4 Then
                CurrentCheck = 23
            ElseIf Form1.Combo1.ListIndex = 5 Then
                CurrentCheck = 24
            End If
            'CurrentCheck = form1.Combo1.ListIndex + 20
            
            If Form1.Combo1.ListIndex = -1 Then
            End If

        End If

        ExeCheckFlag = 1
        If CurrentCheck = 1 Or (XOverList(RelX, RelY).ProgramFlag <> 1 And XOverList(RelX, RelY).ProgramFlag <> 1 + AddNum) Then
            GCSeq1 = Seq1
                    GCSeq2 = Seq2
                    GCSeq3 = Seq3
            
            If XOverList(RelX, RelY).ProgramFlag = 3 And MCTripletFlag = 1 Then
            
            Else
                Call GCCompare
            End If
        End If
        If XOverList(RelX, RelY).ProgramFlag = 0 Or XOverList(RelX, RelY).ProgramFlag = 0 + AddNum Then

            Call RDPChecking
        ElseIf XOverList(RelX, RelY).ProgramFlag = 8 Or XOverList(RelX, RelY).ProgramFlag = 8 + AddNum Then

            Call TSChecking
        ElseIf XOverList(RelX, RelY).ProgramFlag = 1 Or XOverList(RelX, RelY).ProgramFlag = 1 + AddNum Then
            
            
            If pGCTripletflag = 0 Then
                Call GCChecking
            Else
                Call GCChecking2
            End If

        ElseIf XOverList(RelX, RelY).ProgramFlag = 2 Or XOverList(RelX, RelY).ProgramFlag = 2 + AddNum Then

            Call BootscanChecking
        ElseIf XOverList(RelX, RelY).ProgramFlag = 4 Or XOverList(RelX, RelY).ProgramFlag = 4 + AddNum Then

            Call ChimaeraChecking
        ElseIf XOverList(RelX, RelY).ProgramFlag = 5 Or XOverList(RelX, RelY).ProgramFlag = 5 + AddNum Then

            Call SiScanChecking
        ElseIf XOverList(RelX, RelY).ProgramFlag = 6 Or XOverList(RelX, RelY).ProgramFlag = 6 + AddNum Then

            Call PhylProChecking
        ElseIf XOverList(RelX, RelY).ProgramFlag = 3 Or XOverList(RelX, RelY).ProgramFlag = 3 + AddNum Then

            If MCTripletFlag = 0 Then

                Call MaxChiChecking

            Else

                Call MaxChiCheckingII

            End If

        End If

        Dim LineStart As Double
        Dim LineEnd As Double

        If ORFFlag = 1 Then
            Form1.Picture20.BorderStyle = 0
            LineStart = (Form1.Picture20.ScaleWidth / Len(StrainSeq(0))) * RecStart
            LineEnd = (Form1.Picture20.ScaleWidth / Len(StrainSeq(0))) * RecEnd

            If LineEnd > LineStart Then
                Form1.Picture20.DrawMode = 15
                Form1.Picture20.Line (LineStart, 0)-(LineEnd, Form1.Picture20.ScaleHeight), RGB(192, 64, 64), BF
                Form1.Picture20.DrawMode = 13
            Else
                Form1.Picture20.DrawMode = 15
                Form1.Picture20.Line (0, 0)-(LineEnd, Form1.Picture20.ScaleHeight), RGB(192, 64, 64), BF
                Form1.Picture20.Line (LineStart, 0)-(Form1.Picture20.ScaleWidth, Form1.Picture20.ScaleHeight), RGB(192, 64, 64), BF
                Form1.Picture20.DrawMode = 13
            End If

        End If
        
        If LongWindedFlag = 1 Then
            Call UnModNextno
            Call UnModSeqNum(0)
           X = X
            'Seq1 = TreeTraceSeqs(0, Seq1)
            'Seq2 = TreeTraceSeqs(0, Seq2)
            'Seq3 = TreeTraceSeqs(0, Seq3)
            
            'For X = 1 To Len(StrainSeq(0))
            '    SeqNum(X, Seq1) = tSeqNum(X, 0)
            '    SeqNum(X, Seq2) = tSeqNum(X, 1)
            '    SeqNum(X, Seq3) = tSeqNum(X, 2)
            'Next X
        End If
   ' End If
    
    
    Seq1 = oSeq1
    Seq2 = oSeq2
    Seq3 = oSeq3
        
     If XOverList(RelX, RelY).ProgramFlag <> 1 And XOverList(RelX, RelY).ProgramFlag <> 1 + AddNum Then
            GCSeq1 = Seq1
                    GCSeq2 = Seq2
                    GCSeq3 = Seq3
            
            If XOverList(RelX, RelY).ProgramFlag = 3 And MCTripletFlag = 1 Then
            
            Else
                Call GCCompare
            End If
        End If
        
    Form1.Combo1.Enabled = True
    If CurrentCheck <> 20 And CurrentCheck <> 6 And CurrentCheck <> -1 And CurrentCheck <> 21 And CurrentCheck <> 22 And CurrentCheck <> 24 Then
        Form1.Command29.Enabled = True
    End If
    ExeCheckFlag = 0
    Form1.ProgressBar1 = 0
    Form1.SSPanel1.Caption = ""

End Sub
Public Sub OptTree()

Outie = Nextno - 1
Dim LastBoot As Long, TopSeq As Long
    
    LastBoot = 0
    TopBoot = 1
    OldBoot = 0
    cycle = 0
    TSeqLen = Len(StrainSeq(0))
    Do
        For z = 44 To Nextno
            OldBoot = 0
            ReDim ETSeqNum(Len(StrainSeq(0)), Nextno - 1)
            A = -1
            For Y = 0 To Nextno
                If Y <> z Then
                    A = A + 1
                    For X = 1 To Len(StrainSeq(0))
                    
                        ETSeqNum(X, A) = SeqNum(X, Y)
                    Next X
                End If
            Next Y
            Call GetOutie
            Nextno = Nextno - 1
            NumberOfSeqs = Nextno
            Call NJTree2(1)
            Nextno = Nextno + 1
            If TopBoot < OldBoot Then
                TopBoot = OldBoot
                TopSeq = z
            End If
            Form1.SSPanel14.Caption = Str(cycle) + ":" + left(Str(TopBoot / (Nextno + 1)), 6) + ":" + Str(TopSeq) + ":" + Str(z)
        Next z
        If LastBoot < TopBoot Then
            LastBoot = TopBoot
            z = TopSeq
            A = -1
            For Y = 0 To Nextno
                If Y <> z Then
                    A = A + 1
                    StraiName(A) = StraiName(Y)
                    For X = 1 To Len(StrainSeq(0))
                        ETSeqNum(X, A) = SeqNum(X, Y)
                    Next X
                End If
            Next Y
            ReDim SeqNum(Len(StrainSeq(0)), Nextno - 1)
            Nextno = Nextno - 1
            For Y = 0 To Nextno
                'If Y <> Z Then
                    'A = A + 1
                    For X = 1 To Len(StrainSeq(0))
                    
                        SeqNum(X, Y) = ETSeqNum(X, Y)
                    Next X
                'End If
            Next Y
        Else
            Exit Do
        End If
        
        XX = Len(StrainSeq(0))
        ReDim StrainSeq(Nextno)
        For X = 0 To Nextno
            For Y = 1 To XX
                StrainSeq(X) = StrainSeq(X) + Chr(SeqNum(Y, X) - 1)
            Next Y
        Next X
        
        XX = CurDir
        Open "opttree" + Str(cycle) For Output As #1
        Print #1, TopSeq
        For X = 0 To Nextno
            Print #1, ">" + StraiName(X)
            Print #1, StrainSeq(X)
        Next X
        Close #1
        cycle = cycle + 1
        Outie = Nextno - 1
    Loop
    
    XX = Len(StrainSeq(0))
        ReDim StrainSeq(Nextno)
        For X = 0 To Nextno
            For Y = 1 To XX
                StrainSeq(X) = StrainSeq(X) + Chr(SeqNum(Y, X) - 1)
            Next Y
        Next X
        
        XX = CurDir
        Open "opttree" + Str(cycle) For Output As #1
        Print #1, TopSeq
        For X = 0 To Nextno
            Print #1, ">" + StraiName(X)
            Print #1, StrainSeq(X)
        Next X
        Close #1
    
    X = X
End Sub
Public Sub DisableInterface()
    Form1.Combo1.Enabled = False
    Form1.Command29.Enabled = False
    Form1.SSPanel2.Enabled = False
    Form1.Picture8.Enabled = False
    Form1.Picture10.Enabled = False
    Form1.Picture5.Enabled = False
    Form1.Command6.Enabled = False
    Form1.Command25.Enabled = False
End Sub
Public Sub EnableInterface()
    Form1.Combo1.Enabled = True
    Form1.Command29.Enabled = True
    Form1.SSPanel2.Enabled = True
    Form1.Picture8.Enabled = True
    Form1.Picture10.Enabled = True
    Form1.Picture5.Enabled = True
    Form1.Command6.Enabled = True
    Form1.Command25.Enabled = True
End Sub
Public Sub RescanBest(DF)

Dim rCircularflag As Byte, WeightMod() As Long, Scratch() As Integer, SStart As Long, SEnd As Long
'XX = UBound(MaxXOP, 1)
'XX = AddNum
'Exit Sub
rCircularflag = 1
Form1.Timer1.Enabled = False
TManFlag = -1
If DF = 0 Then
    SStart = 1
    SEnd = SEventNumber
    SPF = 1
Else
    SStart = DF
    SEnd = DF
    SPF = 1
End If
oRelX = RelX
oRelY = RelY
For X = SStart To SEnd

   

    Form1.SSPanel1.Caption = "Rescanning" + Str(X) + " of" + Str(SEventNumber) + " detected events"
    RelX = BestEvent(X, 0)
    RelY = BestEvent(X, 1)
    
    If ((XOverList(RelX, RelY).MissIdentifyFlag <> 3 And XOverList(RelX, RelY).MissIdentifyFlag <> 13 And XOverList(RelX, RelY).Accept <> 2) Or (XOverList(RelX, RelY).Accept = 1)) Or DF = 1 Then
        
        xMCTripletFlag = MCTripletFlag
        xMCProportionFlag = MCProportionFlag
        xMCStart = MCStart
        xMCEnd = MCEnd
        xMCMaxP = MCMaxP
        xMCSteplen = MCSteplen
        xMCWinSize = MCWinSize
        xMCWinFract = MCWinFract
        xMCStripGapsFlag = MCStripGapsFlag
        xMCFlag = MCFlag
        xXOverWindowX = XOverWindowX
        xCircularFlag = CircularFlag
        xSpacerFlag = SpacerFlag
        xLowestProb = LowestProb
        xBSTypeFlag = BSTypeFlag
        xBSStepSize = BSStepSize
        xBSStepWin = BSStepWin
        xBSBootReps = BSBootReps
        xBSCutoff = BSCutOff
        xBSPValFlag = BSPValFlag
        xSSGapFlag = SSGapFlag
        xSSVarPFlag = SSVarPFlag
        xSSOutlyerFlag = SSOutlyerFlag
        xSSRndSeed = SSRndSeed
        xSSWinLen = SSWinLen
        xSSStep = SSStep
        xSSNumPerms = SSNumPerms
        xGCMissmatchPen = GCMissmatchPen
        xGCIndelFlag = GCIndelFlag
        xGCMinFragLen = GCMinFragLen
        xGCMinPolyInFrag = GCMinPolyInFrag
        xGCMinPairScore = GCMinPairScore
        xGCMaxOverlapFrags = GCMaxOverlapFrags
        xGCTripletFlag = GCtripletflag
        xPPStripGaps = PPStripGaps
        Dim A As Integer
        If Nextno > UBound(pMaskSeq, 1) Then
            ReDim Preserve pMaskSeq(Nextno)
        End If
        For A = 0 To Nextno
            MaskSeq(A) = pMaskSeq(A)
        Next 'A
        
        MCTripletFlag = pMCTripletFlag
        MCProportionFlag = pMCProportionFlag
        MCStart = pMCStart
        MCEnd = pMCEnd
        MCMaxP = pMCMaxP
        MCSteplen = pMCSteplen
        MCWinSize = pMCWinSize
        MCWinFract = pMCWinFract
        MCStripGapsFlag = pMCStripGapsFlag
        MCFlag = pMCFlag
        XOverWindowX = pXOverWindowX
        CircularFlag = pCircularFlag
        SpacerFlag = pSpacerFlag
        LowestProb = pLowestProb
        BSTypeFlag = pBSTypeFlag
        BSStepSize = pBSStepSize
        BSStepWin = pBSStepWin
        BSBootReps = pBSBootReps
        BSCutOff = pBSCutoff
        BSPValFlag = pBSPValFlag
        
        SSGapFlag = pSSGapFlag
        SSVarPFlag = pSSVarPFlag
        SSOutlyerFlag = pSSOutlyerFlag
        SSRndSeed = pSSRndSeed
        SSWinLen = pSSWinLen
        SSStep = pSSStep
        SSNumPerms = pSSNumPerms
        PPStripGaps = pPPStripGaps
        
        Seq1 = XOverList(RelX, RelY).MajorP
        Seq2 = XOverList(RelX, RelY).MinorP
        Seq3 = XOverList(RelX, RelY).Daughter
        
        ReDim tSeqNum(Len(StrainSeq(0)), 2)
        EN = XOverList(RelX, RelY).Eventnumber
                    
        NJFlag = 0
        Call ModSeqNum(SPF)
                        
                        
        Call MakeTreeSeqs(XOverList(RelX, RelY).Beginning, XOverList(RelX, RelY).Ending)
                        
        GCSeq1 = Seq1
        GCSeq2 = Seq2
        GCSeq3 = Seq3
        If AbortFlag = 1 Then
            
            Exit Sub
        End If
        If rCircularflag = 1 Then
            For z = 0 To Nextno
                MissingData(1, z) = 1
                MissingData(Len(StrainSeq(0)), z) = 1
            Next z
        End If
           
        Call GCCompare
        
        Dim XPS As Long, YPS As Long
        XPS = BestEvent(X, 0)
        YPS = BestEvent(X, 1)
        
        For CurrentCheck = 0 To 10
            GoOn = 200
            If CurrentCheck = 0 And (Confirm(X, 0) = 0 Or XOverList(XPS, YPS).ProgramFlag = 0) Then
                If XOverList(XPS, YPS).ProgramFlag = 0 Then
                    GoOn = 0
                End If
                xSpacerFlag = SpacerFlag
                XOverWindowX = CDbl(Form3.Text2.Text)
                
                Call XOverIII(SPF)
                If SPF = 0 Then Form1.Picture7.Refresh
                SpacerFlag = xSpacerFlag
            ElseIf CurrentCheck = 1 And (Confirm(X, 1) = 0 Or XOverList(XPS, YPS).ProgramFlag = 1) Then
                    If XOverList(XPS, YPS).ProgramFlag = 1 Then
                        GoOn = 1
                    End If
                    Call GCCheck(SPF)
                    If SPF = 0 Then Form1.Picture7.Refresh
                X = X
            
            ElseIf CurrentCheck = 2 And (Confirm(X, 2) = 0 Or XOverList(XPS, YPS).ProgramFlag = 2) Then
                If XOverList(XPS, YPS).ProgramFlag = 2 Then
                    GoOn = 2
                End If
                Call FindSubSeqBS
                
                If DoScans(0, 2) = 1 And UBound(BSFilePos, 2) > 0 Then
        
                    Call BSXoverL(SPF)
                    If SPF = 0 Then Form1.Picture7.Refresh
                Else
                    ReDim Preserve MaxXOP(AddNum - 1, Nextno)
                    
                        '
                        s1col = Yellow
                        s1colb = LYellow
                        s2col = Purple
                        s2colb = LPurple
                        s3col = green
                        s2colb = LGreen
                        Call FindSubSeqBS
                        
                        ReDim Scratch(BSStepWin), WeightMod(BSBootReps, BSStepWin - 1)
                        Dummy = SEQBOOT2(BSRndNumSeed, BSBootReps, BSStepWin, Scratch(0), WeightMod(0, 0))
                        Call BSXoverM(SPF, 0, WeightMod())
                        If SPF = 0 Then Form1.Picture7.Refresh
                        
                    
        
                End If
                
            
                
            ElseIf CurrentCheck = 4 And (Confirm(X, 3) = 0 Or XOverList(XPS, YPS).ProgramFlag = 3) Then
                If XOverList(XPS, YPS).ProgramFlag = 3 Then
                    GoOn = 3
                End If
                Call MCXoverG(SPF)
                If SPF = 0 Then Form1.Picture7.Refresh
                
            ElseIf CurrentCheck = 5 And (Confirm(X, 5) = 0 Or XOverList(XPS, YPS).ProgramFlag = 5) Then
                If XOverList(XPS, YPS).ProgramFlag = 5 Then
                    GoOn = 5
                End If
                Call SSXoverB(SPF)
               If SPF = 0 Then Form1.Picture7.Refresh
               
            ElseIf CurrentCheck = 6 And (Confirm(X, 6) = 0) And (DoScans(0, 7) = 1 Or DoScans(1, 7) = 1) Then
               Call LXoverB(1, 1)
            ElseIf CurrentCheck = 10 And (Confirm(X, 4) = 0 Or XOverList(XPS, YPS).ProgramFlag = 4) Then
                If XOverList(XPS, YPS).ProgramFlag = 4 Then
                    GoOn = 4
                End If
                Call CXoverC(SPF)
                
                If SPF = 0 Then Form1.Picture7.Refresh
            End If
            If GoOn < 200 Then
                Dim RL As Long, LS As Long, BPos As Long, EPos As Long
                If GoOn = 0 Then
                    RL = XOverWindow
                    
                ElseIf GoOn = 1 Then
                    RL = 10
                ElseIf GoOn = 2 Then
                    RL = 10
                ElseIf GoOn = 3 Then
                    RL = MCWinSize / 2
                ElseIf GoOn = 4 Then
                    RL = CWinSize / 2
                ElseIf GoOn = 5 Then
                    RL = 10
                End If
                LS = LenXoverSeq
                BPos = XOverList(XPS, YPS).Beginning
                EPos = XOverList(XPS, YPS).Ending
                
                warn = 0
                Call CheckEndsVB(RL, warn, LS, 0, CircularFlag, Seq1, Seq2, Seq3, BPos, EPos, SeqNum(), XPosDiff(), XDiffpos())
                
                If warn = 1 Then
                    If XOverList(XPS, YPS).SBPFlag = 0 Then
                        XOverList(XPS, YPS).SBPFlag = 1
                    ElseIf XOverList(XPS, YPS).SBPFlag = 2 Then
                        XOverList(XPS, YPS).SBPFlag = 3
                    
                    End If
                'Else
                '    XOverList(XPS, YPS).SBPFlag = 0
                End If
                
                warn = 0
                Call CheckEndsVB(RL, warn, LS, 1, CircularFlag, Seq1, Seq2, Seq3, BPos, EPos, SeqNum(), XPosDiff(), XDiffpos())

                
                If warn = 1 Then
                    If XOverList(XPS, YPS).SBPFlag = 0 Then
                        XOverList(XPS, YPS).SBPFlag = 2
                    ElseIf XOverList(XPS, YPS).SBPFlag = 1 Then
                        XOverList(XPS, YPS).SBPFlag = 3
                    End If
                
                'Else
                '    XOverList(XPS, YPS).SBPFlag = 0
                End If
            End If
        Next CurrentCheck
        
        Call UnModNextno
        Call UnModSeqNum(0)
        'Retrieve variable states saved before analysis began
        MCTripletFlag = xMCTripletFlag
        MCProportionFlag = xMCProportionFlag
        MCStart = xMCStart
        MCEnd = xMCEnd
        MCMaxP = xMCMaxP
        MCSteplen = xMCSteplen
        MCWinSize = xMCWinSize
        MCWinFract = xMCWinFract
        MCStripGapsFlag = xMCStripGapsFlag
        MCFlag = xMCFlag
        XOverWindowX = xXOverWindowX
        CircularFlag = xCircularFlag
        SpacerFlag = xSpacerFlag
        LowestProb = xLowestProb
        SpacerNo = SpacerNo
        BSTypeFlag = xBSTypeFlag
        BSStepSize = xBSStepSize
        BSStepWin = xBSStepWin
        BSBootReps = xBSBootReps
        BSCutOff = xBSCutoff
        BSPValFlag = xBSPValFlag
        SSGapFlag = xSSGapFlag
        SSVarPFlag = xSSVarPFlag
        SSOutlyerFlag = xSSOutlyerFlag
        SSRndSeed = xSSRndSeed
        SSWinLen = xSSWinLen
        SSStep = xSSStep
        SSNumPerms = xSSNumPerms
        GCMissmatchPen = xGCMissmatchPen
        GCIndelFlag = xGCIndelFlag
        GCMinFragLen = xGCMinFragLen
        GCMinPolyInFrag = xGCMinPolyInFrag
        GCMinPairScore = xGCMinPairScore
        GCMaxOverlapFrags = xGCMaxOverlapFrags
        GCtripletflag = xGCTripletFlag
        PPStripGaps = xPPStripGaps
        
    End If
    DoEvents
    Form1.ProgressBar1.Value = (X / SEventNumber) * 100
    Form1.Refresh
    
Next X
RelX = oRelX
RelY = oRelY
Form1.SSPanel1.Caption = ""
Form1.ProgressBar1.Value = 0
Screen.MousePointer = 0
Form1.Timer1.Enabled = True
End Sub

Public Sub UpdateEvents(SE)

Dim BestP, AProg() As Byte
ReDim AProg(AddNum * 2)

If DoScans(0, 0) = 1 Then AProg(0) = 1: AProg(0 + AddNum) = 1
If DoScans(0, 1) = 1 Then AProg(1) = 1: AProg(1 + AddNum) = 1
If DoScans(0, 2) = 1 Then AProg(2) = 1: AProg(2 + AddNum) = 1
If DoScans(0, 3) = 1 Then AProg(3) = 1: AProg(3 + AddNum) = 1
If DoScans(0, 4) = 1 Then AProg(4) = 1: AProg(4 + AddNum) = 1
If DoScans(0, 5) = 1 Then AProg(5) = 1: AProg(5 + AddNum) = 1
'XX = UBound(BestP, 1)
'ReDim BestP(SEventNumber)
'For Y = 0 To SEventNumber
    For X = 0 To AddNum - 1
        Confirm(SE, X) = 0
        ConfirmP(SE, X) = 0
    Next X
'Next Y
BestP = 100000
    For X = 0 To Nextno
        For Y = 1 To CurrentXOver(X)
            If SuperEventList(XOverList(X, Y).Eventnumber) = SE Then
                pf = XOverList(X, Y).ProgramFlag
                If pf <= AddNum - 1 Then
                    If XOverList(X, Y).Probability > 0 And (XOverList(X, Y).Probability < BestP Or BestP = 0) Then
                        If AProg(XOverList(X, Y).ProgramFlag) = 1 Then
                            BestP = XOverList(X, Y).Probability
                            BestEvent(SE, 0) = X
                            BestEvent(SE, 1) = Y
                        End If
                    End If
                    Confirm(SE, pf) = Confirm(SE, pf) + 1
                    ConfirmP(SE, pf) = ConfirmP(SE, pf) + -Log10(XOverList(X, Y).Probability)
                End If
            End If
        Next Y
    Next X
    If BestP = 100000 Then
        For X = 0 To Nextno
            For Y = 1 To CurrentXOver(X)
                If SuperEventList(XOverList(X, Y).Eventnumber) = SE Then
                    pf = XOverList(X, Y).ProgramFlag
                    If pf > AddNum - 1 Then
                        If XOverList(X, Y).Probability > 0 And (XOverList(X, Y).Probability < BestP Or BestP = 0) Then
                            If AProg(XOverList(X, Y).ProgramFlag) = 1 Then
                                BestP = XOverList(X, Y).Probability
                                BestEvent(SE, 0) = X
                                BestEvent(SE, 1) = Y
                            End If
                        End If
                    
                        pf = pf - AddNum
                        Confirm(SE, pf) = Confirm(SuperEventList(XOverList(X, Y).Eventnumber), pf) + 1
                        ConfirmP(SE, pf) = ConfirmP(SuperEventList(XOverList(X, Y).Eventnumber), pf) + -Log10(XOverList(X, Y).Probability)
                    End If
                End If
            Next Y
        Next X
    End If
End Sub

Public Sub CheckDists(SeqNum() As Integer, PermValid() As Double)

Dim PVA() As Long
    ReDim PVA(Nextno, Nextno)
    For X = 0 To Nextno
        For Y = X + 1 To Nextno
            For z = 1 To Len(StrainSeq(0))
                If SeqNum(z, X) <> 46 And SeqNum(z, Y) <> 46 Then
                    PVA(X, Y) = PVA(X, Y) + 1
                End If
            Next z
            PVA(Y, X) = PVA(X, Y)
        Next Y
        
    Next X
    For X = 0 To Nextno
        For Y = X + 1 To Nextno
            If PVA(X, Y) <> PermValid(X, Y) Then
                X = X
            End If
        Next Y
    Next X
End Sub

Public Sub FindBestRecSignalVB(TotalNoRecombinants, LowP As Double, Trace() As Long, PCurrentXOver() As Integer, DoneSeq() As Byte, PXOList() As XOverDefine)
pvalsearch = 0
If pvalsearch = 1 Then
        For X = 0 To Nextno
            TotalNoRecombinants = TotalNoRecombinants + PCurrentXOver(X)
           
            For Y = 1 To PCurrentXOver(X)
                
             '   If SEventNumber > 68 And X = X Then
             '
             '       If PXOList(X, Y).Daughter = 16 Or PXOList(X, Y).Daughter = 29 Or PXOList(X, Y).Daughter = 15 Then
             '       If PXOList(X, Y).MinorP = 16 Or PXOList(X, Y).MinorP = 29 Or PXOList(X, Y).MinorP = 15 Then
              ''      If PXOList(X, Y).MajorP = 16 Or PXOList(X, Y).MajorP = 29 Or PXOList(X, Y).MajorP = 15 Then
             '           X = X
             '
             '       End If
             '       End If
             '       End If
             '   End If
                If DoneSeq(X, Y) = DoneTarget Then
                    CPVal = PXOList(X, Y).Probability
                    
                    If CPVal > 0 And CPVal < LowP Then ' And (PXOList(X, Y).ProgramFlag <> 5 Or OnlySiScan = 1) Then
                        If PXOList(X, Y).Beginning <> PXOList(X, Y).Ending Then
                            LowP = CPVal
                            Trace(0) = X
                            Trace(1) = Y
                            
                       End If
                    
                       ' XX = PXOList(Trace(0), Trace(1)).Daughter
                       ' XX = PXOList(Trace(0), Trace(1)).MinorP
                       ' XX = PXOList(Trace(0), Trace(1)).MajorP
                       ' XX = PXOList(Trace(0), Trace(1)).Beginning
                       ' XX = PXOList(Trace(0), Trace(1)).Ending
                    End If
                End If
            Next Y
        Next X
    Else 'use a distance search
        Dim DistA() As Double, DistB() As Double, Sites() As Double, LOwdist As Double, s(2) As Long, SubDiffs() As Double, SubValid() As Double
        LOwdist = 1000
        'LowP = 100000
        ReDim DistA(Nextno, Nextno), DistB(Nextno, Nextno), SubValid(Nextno, Nextno), SubDiffs(Nextno, Nextno)
        For X = 0 To Nextno
            TotalNoRecombinants = TotalNoRecombinants + PCurrentXOver(X)
           
            For Y = 1 To PCurrentXOver(X)
                
                If DoneSeq(X, Y) = DoneTarget Then
                    s(0) = PXOList(X, Y).Daughter
                    s(1) = PXOList(X, Y).MajorP
                    s(2) = PXOList(X, Y).MinorP
                    BE = PXOList(X, Y).Beginning
                    EN = PXOList(X, Y).Ending
                    
                    
                    Dummy = vQuickDist(Len(StrainSeq(0)), Nextno, Nextno, BE, EN, DistA(0, 0), DistB(0, 0), SubValid(0, 0), SubDiffs(0, 0), PermValid(0, 0), PermDiffs(0, 0), SeqNum(0, 0), s(0))
                    
                    
                    
                    
                    For B = 0 To 2
                        For C = B + 1 To 2
                            If DistA(s(B), s(C)) < LOwdist Then
                                LOwdist = DistA(s(B), s(C))
                                Trace(0) = X
                                Trace(1) = Y
                                LowP = PXOList(X, Y).Probability
                            ElseIf DistA(s(B), s(C)) = LOwdist Then
                                If PXOList(X, Y).Probability < LowP Then
                                    LOwdist = DistA(s(B), s(C))
                                    Trace(0) = X
                                    Trace(1) = Y
                                    LowP = PXOList(X, Y).Probability
                                End If
                            End If
                            If DistB(s(B), s(C)) < LOwdist Then
                                LOwdist = DistB(s(B), s(C))
                                Trace(0) = X
                                Trace(1) = Y
                                LowP = PXOList(X, Y).Probability
                            ElseIf DistB(s(B), s(C)) = LOwdist Then
                                If PXOList(X, Y).Probability < LowP Then
                                    LOwdist = DistB(s(B), s(C))
                                    Trace(0) = X
                                    Trace(1) = Y
                                    LowP = PXOList(X, Y).Probability
                                End If
                            End If
                        Next C
                    Next B
                   
                    
                End If
            Next Y
        Next X
        'LowP = PXOList(Trace(0), Trace(1)).Probability
    End If
End Sub
'Declare Function FindActualEvents Lib "dna.dll" (RSize As Long, ByRef Don As Byte, ByRef BPMatch As Long,  ByRef OKSeq As Double, ByRef FoundOne As Long, ByRef SP As Long, ByRef EP As Long, ByRef RCorr As Double, ByRef OLSeq As Long, ByRef OLSeqB As Long, ByRef OLSeqE As Long, ByRef CSeq As LoInvS As Byte, ByRef TrS As Long, ByRef TMatch As Double, ByRef PXOList As XOverDefine, ByRef PCurrentXOver As Integer, ByRef SQ As Long, ByRef tDon As Byte, ByRef ISeqs As Long, ByRef TList As Byte, ByRef CompMat As Long) As Long

Public Sub FindActualEventsVB(RLScore() As Double, UNF() As Byte, InvList() As Long, Nextno As Long, RSize() As Long, BPMatch() As Long, BMatch() As Double, OKSeq() As Double, FoundOne() As Long, SP() As Long, EP() As Long, RCorr() As Double, OLSeq() As Long, OLSeqB() As Long, OLSeqE() As Long, CSeq() As Long, RNum() As Long, RList() As Long, InvS() As Byte, TMatch() As Double, PXOList() As XOverDefine, PCurrentXOver() As Integer, SQ() As Long, tDon() As Byte, ISeqs() As Long, CompMat() As Long)
    Dim TList() As Byte, Don() As Byte, TrS(2) As Long, Y As Long, X As Long
        
        
  '      If X = 12345 Then ' use this if all must be compared against all
  '          ReDim TList(2, NextNo)
  '          For X = 0 To 2
  '              For Y = 0 To RNum(X)
  '                  TList(X, RList(X, Y)) = 1
  '              Next Y
  '          Next X
  '      End If
        
        'If SEventNumber = 24 Then
        '    For X = 0 To NextNo
        '        For Y = 1 To PCurrentXOver(X)
        '            If PXOList(X, Y).Daughter = ISeqs(0) Or PXOList(X, Y).MinorP = ISeqs(0) Or PXOList(X, Y).MajorP = ISeqs(0) Then
        '            If PXOList(X, Y).Daughter = ISeqs(1) Or PXOList(X, Y).MinorP = ISeqs(1) Or PXOList(X, Y).MajorP = ISeqs(1) Then
        '            'If PXOList(X, Y).MajorP = ISeqs(0) Or PXOList(X, Y).MajorP = ISeqs(1) Or PXOList(X, Y).MajorP = 13 Then
        '                X = X
        '            'End If
        '            End If
        '            End If
        '        Next Y
        '    Next X
        '    X = X
        'End If
        'X = X
        'For X = 0 To 2
        '    XX = RNum(X) '23,15,15:15,15,23
        '    XX = InPen(X) '1(dug),1(10cd),0(cin)
        '    For Y = 0 To RNum(X)
        '
        '        XX = StraiName(RList(X, Y))
        '        'dug,dug
        '        '10cd
        '        'cin,08bc,ges,cza,cbr
        '    Next Y
        'Next X
    'End If
        
        ReDim Don(2, Nextno)
        
        ReDim tDon(5)
        'XX = Trace(1)
        For WinPP = 0 To 2
            
            'If X = X Then 'use this to only screen against iseqs
                ReDim TList(2, Nextno)
                TList(CompMat(WinPP, 0), ISeqs(CompMat(WinPP, 0))) = 1
                TList(CompMat(WinPP, 1), ISeqs(CompMat(WinPP, 1))) = 1
                For Y = 0 To RNum(WinPP)
                    TList(WinPP, RList(WinPP, Y)) = 1
                Next Y
            'End If
            
            ReDim FoundOne(Nextno)
            If X = X Then
                Dummy = FindActualEvents(Len(StrainSeq(0)), WinPP, Nextno, UBound(PXOList, 1), RSize(0), Don(0, 0), BPMatch(0, 0, 0), BMatch(0, 0), OKSeq(0, 0, 0), FoundOne(0), SP(0), EP(0), RCorr(0, 0, 0), OLSeq(0), OLSeqB(0), OLSeqE(0), CSeq(0), RNum(0), RList(0, 0), InvS(0, 0), TrS(0), TMatch(0), PXOList(0, 0), PCurrentXOver(0), SQ(0), tDon(0), ISeqs(0), TList(0, 0), CompMat(0, 0))
                
            Else
                For X = 0 To Nextno
                    If TList(WinPP, X) = 1 Or ISeqs(CompMat(WinPP, 0)) = X Or ISeqs(CompMat(WinPP, 1)) = X Then
                        OldY = -1
                        For Y = 1 To PCurrentXOver(X)
                            If OldY <> Y Then
                                tDon(0) = 0: tDon(1) = 0: tDon(2) = 0: tDon(3) = 0: tDon(4) = 0: tDon(5) = 0
                                OldY = Y
                                zzx = 0
                            Else
                                zzx = zzx + 1
                                If zzx > 6 Then
                                    ReDim tDon(5)
                                    OldY = Y
                                    Y = Y + 1
                                    zzx = 0
                                    If Y > PCurrentXOver(X) Then Exit For
                                    
                                End If
                                If Y = OldY Then
                                    ZZZX = 0
                                    For A = 0 To 5
                                        ZZZX = ZZZX + tDon(A)
                                    Next A
                                    If ZZZX = 6 Then
                                        Y = Y + 1
                                        ReDim tDon(5)
                                        OldY = Y
                                        zzx = 0
                                        ZZZX = 0
                                        If Y > PCurrentXOver(X) Then Exit For
                                    End If
                                    
                                End If
                            End If
                            SQ(1) = PXOList(X, Y).MajorP
                            If TList(WinPP, SQ(1)) = 1 Or ISeqs(CompMat(WinPP, 0)) = SQ(1) Or ISeqs(CompMat(WinPP, 1)) = SQ(1) Then
                                GoOn = 0
                                SQ(2) = PXOList(X, Y).MinorP
                                If X = 1234 Then
                                        
                                   ' Dummy = tester1(X, Y, GoOn, Len(StrainSeq(0)), WinPP, NextNo, UBound(PXOList, 1), RSize(0), Don(0, 0), BPMatch(0, 0, 0), BMatch(0, 0), OKSeq(0, 0, 0), FoundOne(0), SP(0), EP(0), RCorr(0, 0, 0), OLSeq(0), OLSeqB(0), OLSeqE(0), CSeq(0), RNum(0), RList(0, 0), InvS(0, 0), TrS(0), TMatch(0), PXOList(0, 0), PCurrentXOver(0), SQ(0), tDon(0), ISeqs(0), TList(0, 0), CompMat(0, 0))
                                Else
                                    If (TList(WinPP, SQ(2)) = 1 Or ISeqs(CompMat(WinPP, 0)) = SQ(2) Or ISeqs(CompMat(WinPP, 1)) = SQ(2)) Then
                                        SQ(0) = X
                                     '   If SEventNumber = 24 Then
                                     '   'XX = RNum(1)
                                     '   'XX = StraiName(ISeqs(WinPP))
                                     '       If SQ(0) = ISeqs(0) Or SQ(0) = ISeqs(1) Or SQ(0) = 33 Then
                                     '       If SQ(1) = ISeqs(0) Or SQ(1) = ISeqs(1) Or SQ(1) = 33 Then
                                     '       If SQ(2) = ISeqs(0) Or SQ(2) = ISeqs(1) Or SQ(2) = 33 Then
                                     '           X = X
                                     '       End If
                                     '       End If
                                     '       End If
                                     '   End If
                                        TMatch(0) = 0
                                        If tDon(0) = 0 And ((WinPP = 0 And Don(0, SQ(0)) = 0) Or (WinPP = 1 And Don(1, SQ(1)) = 0) Or (WinPP = 2 And Don(2, SQ(2)) = 0)) And TList(0, SQ(0)) = 1 And TList(1, SQ(1)) = 1 And TList(2, SQ(2)) = 1 Then
                                            tDon(0) = 1: TrS(0) = SQ(0): TrS(1) = SQ(1): TrS(2) = SQ(2): TMatch(0) = 3
                                        ElseIf tDon(1) = 0 And ((WinPP = 0 And Don(0, SQ(0)) = 0) Or (WinPP = 1 And Don(1, SQ(2)) = 0) Or (WinPP = 2 And Don(2, SQ(1)) = 0)) And TList(0, SQ(0)) = 1 And TList(1, SQ(2)) = 1 And TList(2, SQ(1)) = 1 Then
                                            tDon(1) = 1: TrS(0) = SQ(0): TrS(1) = SQ(2): TrS(2) = SQ(1): TMatch(0) = 3
                                        ElseIf tDon(2) = 0 And ((WinPP = 0 And Don(0, SQ(1)) = 0) Or (WinPP = 1 And Don(1, SQ(2)) = 0) Or (WinPP = 2 And Don(2, SQ(0)) = 0)) And TList(0, SQ(1)) = 1 And TList(1, SQ(2)) = 1 And TList(2, SQ(0)) = 1 Then
                                            tDon(2) = 1: TrS(0) = SQ(1): TrS(1) = SQ(2): TrS(2) = SQ(0): TMatch(0) = 3
                                        ElseIf tDon(3) = 0 And ((WinPP = 0 And Don(0, SQ(1)) = 0) Or (WinPP = 1 And Don(1, SQ(0)) = 0) Or (WinPP = 2 And Don(2, SQ(2)) = 0)) And TList(0, SQ(1)) = 1 And TList(1, SQ(0)) = 1 And TList(2, SQ(2)) = 1 Then
                                            tDon(3) = 1: TrS(0) = SQ(1): TrS(1) = SQ(0): TrS(2) = SQ(2): TMatch(0) = 3
                                        ElseIf tDon(4) = 0 And ((WinPP = 0 And Don(0, SQ(2)) = 0) Or (WinPP = 1 And Don(1, SQ(1)) = 0) Or (WinPP = 2 And Don(2, SQ(0)) = 0)) And TList(0, SQ(2)) = 1 And TList(1, SQ(1)) = 1 And TList(2, SQ(0)) = 1 Then
                                            tDon(4) = 1: TrS(0) = SQ(2): TrS(1) = SQ(1): TrS(2) = SQ(0): TMatch(0) = 3
                                        ElseIf tDon(5) = 0 And ((WinPP = 0 And Don(0, SQ(2)) = 0) Or (WinPP = 1 And Don(1, SQ(0)) = 0) Or (WinPP = 2 And Don(2, SQ(1)) = 0)) And TList(0, SQ(2)) = 1 And TList(1, SQ(0)) = 1 And TList(2, SQ(1)) = 1 Then
                                            tDon(5) = 1: TrS(0) = SQ(2): TrS(1) = SQ(0): TrS(2) = SQ(1): TMatch(0) = 3
                                        End If
                                        
                                        If TMatch(0) = 3 Then
                                        
                                            For A = 0 To 2
                                                If TrS(WinPP) = SQ(A) Then
                                                    
                                                    Exit For
                                                End If
                                            Next A
                                        
                                        
                                            If InvS(CompMat(WinPP, 0), TrS(CompMat(WinPP, 0))) = 0 And InvS(CompMat(WinPP, 1), TrS(CompMat(WinPP, 1))) = 0 Then
                                                For A = 0 To RNum(WinPP)
                                                    If RList(WinPP, A) = TrS(WinPP) Then
                                                        Exit For
                                                    End If
                                                Next A
                                                
                                                If A > RNum(WinPP) Then
                                                    TMatch(0) = 0
                                                Else
                                                    CSeq(1) = A
                                                    GoOn = 1
                                                    'If FoundOne(CSeq(1)) = 1 Then TMatch(0) = 0
                                                End If
                                            Else
                                                TMatch(0) = 0
                                            End If
                                            
                                            'ie an event involving a potentially recombinant sequence is found.
                                                            
                                                If TMatch(0) = 3 And GoOn = 1 Then
                                                    'check for region overlap
                                                    
                                                    BPos2 = PXOList(X, Y).Beginning
                                                    EPos2 = PXOList(X, Y).Ending
                                                    OLSize = FindOverlap(Len(StrainSeq(0)), BPos2, EPos2, RSize(0), OLSeq(0))
                                                    If OLSize > 0 Then
                                                        TMatch(1) = (OLSize * 2) / (RSize(0) + RSize(1))
                                                    Else
                                                        TMatch(1) = 0
                                                    End If
                                                    otMatch = TMatch(1)
                                                    If TMatch(0) * TMatch(1) > 1 Then
                                                       
                                                            'its in the right region but is it the same event?
                                                            If RCorr(WinPP, 2, RList(WinPP, CSeq(1))) > 0.83 And TMatch(1) > 0.6 Then
                                                                TMatch(0) = 1
                                                            End If
                                                            
                                                            If RCorr(WinPP, 2, RList(WinPP, CSeq(1))) > 0.83 Or RCorr(WinPP, 0, RList(WinPP, CSeq(1))) > 0.83 Then
                                                                OLSize = FindOverlap(Len(StrainSeq(0)), SP(0), EP(1), RSize(2), OLSeqB(0))
                                                                If OLSize > 0 Then
                                                                    TMatch(1) = (OLSize * 2) / (RSize(2) + RSize(3))
                                                                Else
                                                                    TMatch(1) = 0
                                                                End If
                                                                If TMatch(1) > 0.2 Then
                                                                    TMatch(0) = TMatch(0) + 1
                                                                ElseIf RCorr(WinPP, 0, RList(WinPP, CSeq(1))) > 0.83 Then
                                                                    If TMatch(1) = 0 Or OLSize = RSize(2) Then
                                                                        TMatch(0) = TMatch(0) - 0.5
                                                                    End If
                                                                End If
                                                                
                                                            End If
                                                            
                                                            If RCorr(WinPP, 2, RList(WinPP, CSeq(1))) > 0.83 Or RCorr(WinPP, 1, RList(WinPP, CSeq(1))) > 0.83 Then
                                                                OLSize = FindOverlap(Len(StrainSeq(0)), SP(2), EP(3), RSize(4), OLSeqE(0))
                                                                If OLSize > 0 Then
                                                                    TMatch(1) = (OLSize * 2) / (RSize(4) + RSize(5))
                                                                Else
                                                                    TMatch(1) = 0
                                                                End If
                                                                
                                                                If TMatch(1) > 0.2 Then
                                                                    TMatch(0) = TMatch(0) + 1
                                                                ElseIf RCorr(WinPP, 0, RList(WinPP, CSeq(1))) > 0.83 Then
                                                                    If TMatch(1) = 0 Or OLSize = RSize(4) Then
                                                                        TMatch(0) = TMatch(0) - 0.5
                                                                    End If
                                                                End If
                                                                
                                                            End If
                                                            'swap around parents and recombinants later (after the next detection step)
                                                            If TMatch(0) >= 1 Then
                                                                FoundOne(CSeq(1)) = 1
                                                                'If RList(0, CSeq(1)) = 20 Then
                                                                '    X = X
                                                                'End If
                                                                If BMatch(WinPP, RList(WinPP, CSeq(1))) < otMatch Then
                                                                    OKSeq(WinPP, 1, RList(WinPP, CSeq(1))) = otMatch
                                                                    BMatch(WinPP, RList(WinPP, CSeq(1))) = otMatch
                                                                    BPMatch(WinPP, 0, RList(WinPP, CSeq(1))) = BPos2
                                                                    BPMatch(WinPP, 1, RList(WinPP, CSeq(1))) = EPos2
                                                                End If
                                                                Y = Y - 1
                                                            End If
                                                        
                                                        
                                                       
                                                    End If
                                                Else
                                                    Y = Y - 1
                                                End If
                                        End If
                                    End If
                                End If
                            End If
                        Next Y
                    End If
                Next X
            End If
            For X = 0 To RNum(WinPP)
                If FoundOne(X) = 0 Or InvList(WinPP, X) = 1 Then
                    UNF(WinPP, RList(WinPP, X)) = 1
                End If
            Next X
            
            Dummy = StripUnfound(WinPP, RNum(0), RList(0, 0), InvList(0, 0), FoundOne(0), RCorr(0, 0, 0), RLScore(0, 0))
           
        Next WinPP
End Sub

Public Sub SignalCount(XOverList() As XOverDefine, CurrentXOver() As Integer)

Dim TotE() As Long, Pr As Integer, Prob As Double
ReDim TotE(30), oRecombNo(100)
    
For X = 0 To Nextno
    For Y = 1 To CurrentXOver(X)
        Pr = XOverList(X, Y).ProgramFlag
        Prob = XOverList(X, Y).Probability
        If Prob < LowestProb And Prob > 0 Then
            oRecombNo(Pr) = oRecombNo(Pr) + 1
        End If
    Next Y
Next X
For X = 0 To AddNum - 1
    oRecombNo(100) = oRecombNo(100) + oRecombNo(X)
Next X

End Sub
Public Sub Finalise(BCurrentXOver() As Integer, BestXOList() As XOverDefine)
Dim SO As String, MXOSize As Long, X As Long, Y As Long

Nextno = oNextno
AbortFlag = 0
MXOSize = 0
For X = 0 To Nextno
    If MXOSize < BCurrentXOver(X) Then MXOSize = BCurrentXOver(X)
Next X

XOSize = MXOSize + 10

ReDim CurrentXOver(Nextno), XOverList(Nextno, XOSize)

For X = 0 To Nextno
    CurrentXOver(X) = BCurrentXOver(X)
    If UBound(BestXOList, 2) < CurrentXOver(X) Then
        ReDim Preserve BestXOList(UBound(BestXOList, 1), CurrentXOver(X))
    End If
    For Y = 1 To CurrentXOver(X)
        XOverList(X, Y) = BestXOList(X, Y)
    Next Y
Next X

Call SignalCount(XOverList(), CurrentXOver())
Call UpdateRecNums(SEventNumber)

ReDim SeqNum(Len(StrainSeq(0)), Nextno + 1)
For X = 0 To Nextno
    StraiName(X) = PermStraiName(X)
Next X

Dummy = CopySeqs(Len(StrainSeq(0)), UBound(PermSeqNum, 2), PermSeqNum(0, 0), SeqNum(0, 0))


ReDim TreeDistance(Nextno, Nextno)
ReDim Preserve SeqCol(Nextno), FFillCol(Nextno)
ReDim Preserve XCord(4, 3, Nextno + 2), YCord(4, 3, Nextno + 2), RYCord(4, 3, Nextno + 2)
ReDim Distance(Nextno, Nextno), PermValid(Nextno, Nextno), PermDiffs(Nextno, Nextno)
ReDim Preserve StraiName(Nextno)

UDst = PermUDst
AvDst = PermAvDst
For X = 0 To PermNextNo
    For Y = 0 To PermNextNo
        Distance(X, Y) = PermDistance(X, Y)
        TreeDistance(X, Y) = PermTreeDistance(X, Y)
        PermValid(X, Y) = PPermValid(X, Y)
        PermDiffs(X, Y) = PPermDiffs(X, Y)
    Next Y
Next X
LSeq = Len(StrainSeq(0))

TreeDistFlag = 1
NJFlag = 1
'GCIndelFlag = oGCIndelFlag
Call UPGMA(0, 0)
Call DrawTree

Call SetUpEvents

oRecombNo(100) = 0
UpdateRecNums (SEventNumber)
ET = Abs(GetTickCount)
Form1.Label57(0).Caption = DoTimeII(Abs(ET - ST))
oRecombNo(100) = SEventNumber
End Sub

Public Sub Rewind(TestC() As Long, BCurrentXOver() As Integer, BestXOList() As XOverDefine)

Call Finalise(BCurrentXOver(), BestXOList())
For X = 0 To Nextno
    For Y = 1 To CurrentXOver(X)
        If TestC(SuperEventList(XOverList(X, Y).Eventnumber)) <> 1 Or SuperEventList(XOverList(X, Y).Eventnumber) = SEventNumber Then
            If XOverList(X, Y).Accept <> 1 Then
                XOverList(X, Y).Accept = 3
            End If
        End If
    Next Y
Next X
For X = 0 To Nextno
    For Y = 1 To BCurrentXoverMi(X)
        If TestC(SuperEventList(BestXOListMi(X, Y).Eventnumber)) <> 1 Or SuperEventList(BestXOListMi(X, Y).Eventnumber) = SEventNumber Then
            If BestXOListMi(X, Y).Accept <> 1 Then
                BestXOListMi(X, Y).Accept = 3
            End If
        End If
    Next Y
Next X
For X = 0 To Nextno
    For Y = 1 To BCurrentXoverMa(X)
        If TestC(SuperEventList(BestXOListMa(X, Y).Eventnumber)) <> 1 Or SuperEventList(BestXOListMa(X, Y).Eventnumber) = SEventNumber Then
            If BestXOListMa(X, Y).Accept <> 1 Then
                BestXOListMa(X, Y).Accept = 3
            End If
        End If
    Next Y
Next X
End Sub
Public Sub RemoveAccepts(BCurrentXOver() As Integer, BestXOList() As XOverDefine)

For X = 0 To PermNextNo
    For Y = 1 To BCurrentXOver(X)
        If BestXOList(X, Y).Accept = 3 Then
            If BestXOList(X, Y).ProgramFlag < AddNum Then
                BestXOList(X, Y).Accept = 0
            Else
                BestXOList(X, Y).Accept = 2
            End If
        End If
    Next Y
Next X
For X = 0 To PermNextNo
    For Y = 1 To BCurrentXoverMi(X)
        If BestXOListMi(X, Y).Accept = 3 Then
            If BestXOListMi(X, Y).ProgramFlag < AddNum Then
                BestXOListMi(X, Y).Accept = 0
            Else
                BestXOListMi(X, Y).Accept = 2
            End If
        End If
    Next Y
Next X
For X = 0 To PermNextNo
    For Y = 1 To BCurrentXoverMa(X)
        If BestXOListMa(X, Y).ProgramFlag < AddNum Then
                BestXOListMa(X, Y).Accept = 0
            Else
                BestXOListMa(X, Y).Accept = 2
            End If
    Next Y
Next X
End Sub

Public Sub TestConflict(TestC() As Long, WinPP, CurAge As Double, Fail, BPos, EPos, TraceSub() As Long, AgeEvent() As Double, AgeScore() As Double, EventScore() As Long, ISeqs() As Long, BestXOList() As XOverDefine, BCurrentXOver() As Integer)
Dim X As Long, Y As Long, z As Long, ProbPar As Long, EarliestEvent As Long
'curage = ageevent(1)
Fail = 0
EarliestEvent = SEventNumber

ReDim TestC(SEventNumber)
Dim TestB()
ReDim TestB(SEventNumber)
For X = 0 To PermNextNo
    For Y = 1 To BCurrentXOver(X)
        XX = SuperEventList(BestXOList(X, Y).Eventnumber)
        
        If BestXOList(X, Y).DHolder < 0 Then
            ''If XX = 7 Then
            '    X = X
            'End If
            'If EarliestEvent > SuperEventList(BestXOList(X, Y).Eventnumber) Then
                EarliestEvent = SuperEventList(BestXOList(X, Y).Eventnumber)
                If BestXOList(X, Y).MinorP = ISeqs(WinPP) Or BestXOList(X, Y).MajorP = ISeqs(WinPP) Then
                    If BestXOList(X, Y).MinorP = ISeqs(WinPP) Then
                        ProbPar = BestXOList(X, Y).MinorP
                    Else
                        ProbPar = BestXOList(X, Y).MajorP
                    End If
                    BPos2 = BestXOList(X, Y).Beginning
                    EPos2 = BestXOList(X, Y).Ending
                    If BPos2 < EPos2 Then
                        For z = BPos2 To EPos2
                            If AgeScore(z, ProbPar) > -1 Then
                                ' AgeEvent(1, SEventNumber)
                                If AgeScore(z, ProbPar) <= CurAge Then
                                    TestC(EarliestEvent) = TestC(EarliestEvent) + 1
                                Else
                                    TestB(EarliestEvent) = TestB(EarliestEvent) + 1
                                End If
                            End If
                        Next z
                    Else
                        For z = 1 To EPos2
                            If AgeScore(z, ProbPar) > -1 Then
                                If AgeScore(z, ProbPar) <= CurAge Then
                                     TestC(EarliestEvent) = TestC(EarliestEvent) + 1
                                Else
                                    TestB(EarliestEvent) = TestB(EarliestEvent) + 1
                                End If
                            End If
                        Next z
                        For z = BPos2 To Len(StrainSeq(0))
                            If AgeScore(z, ProbPar) > -1 Then
                                If AgeScore(z, ProbPar) < CurAge Then
                                     TestC(EarliestEvent) = TestC(EarliestEvent) + 1
                                Else
                                    TestB(EarliestEvent) = TestB(EarliestEvent) + 1
                                End If
                            End If
                        Next z
                    End If
                End If
            'End If
        End If
    Next Y
Next X


'For X = 0 To PermNextNo
'    For Y = 1 To BCurrentXOver(X)
'        If BestXOList(X, Y).DHolder < 0 Then
'            If SuperEventList(BestXOList(X, Y).Eventnumber) = 7 Then
'                XX = BestXOList(X, Y).Beginning '58-1065
'                XX = BestXOList(X, Y).Ending
'                XX = BPos
'                XX = EPos
'            End If
'        End If
'    Next Y
'Next X

'maybe handle events that overlap the edges of recombinant regions separately
'44-1512,2066-2182,1204-1262,1263-1281,1498-1512,1377-1985
If BPos < EPos Then
    MinR = (EPos - BPos) / 2
Else
    MinR = (Len(StrainSeq(0)) - BPos + EPos) / 2
End If
MinR = 0 'any overlap of a newer event is not allowed
EarliestEvent = SEventNumber + 1
GoOn = 0
For X = 1 To SEventNumber
    If TestC(X) > 0 Then
        If TestC(X) > MinR Or (TestC(X) / (TestB(X) + TestC(X))) > 0.5 Then
            TestC(X) = 1
            GoOn = 1
            If X < EarliestEvent Then
                EarliestEvent = X
            End If
        Else
            TestC(X) = 0
        End If
    End If

Next X

If GoOn = 0 Then Exit Sub
'For X = EarliestEvent To SEventNumber - 1
'    TestC(X) = 1
'Next X

Dim ExtremeTrim As Byte, DE() As Long, CP As Long, BPos3 As Long, EPos3 As Long, FragReg() As Byte, CTE As Long
ReDim DE(SEventNumber)
DE(SEventNumber) = 1
Dim EE As Long
EE = EarliestEvent
ExtremeTrim = 1 'remove all later events with reference to any sequence in which an event has been removed.
Dim YH As Long



For X = 0 To PermNextNo
    For Y = 1 To BCurrentXOver(X)
        
        If BestXOList(X, Y).DHolder < 0 Then
            CTE = SuperEventList(BestXOList(X, Y).Eventnumber)
            If CTE >= EE Then '5
                If DE(CTE) = 0 Then
                    If TestC(CTE) = 1 Then
                        YH = BCurrentXOver(X)
                        'DE(CTE) = 1
                        CP = BestXOList(X, Y).Daughter
                        
                        If ExtremeTrim = 0 Then
                            BPos3 = BestXOList(X, Y).Beginning
                            EPos3 = BestXOList(X, Y).Ending
                            ReDim FragReg(Len(StrainSeq(0)))
                            If BPos3 < EPos3 Then
                                For z = BPos3 To EPos3
                                    FragReg(z) = 1
                                Next z
                            Else
                                For z = 1 To EPos3
                                    FragReg(z) = 1
                                Next z
                                For z = BPos3 To Len(StrainSeq(0))
                                    FragReg(z) = 1
                                Next z
                            End If
                        End If
                        CurAge = AgeEvent(1, CTE)
                        For A = 0 To PermNextNo
                            For B = 1 To BCurrentXOver(A)
                                If DE(SuperEventList(BestXOList(A, B).Eventnumber)) = 0 Then
                                    'If BestXOList(A, B).DHolder < 0 Then
                                        If SuperEventList(BestXOList(A, B).Eventnumber) > CTE Then
                                            If TestC(SuperEventList(BestXOList(A, B).Eventnumber)) = 0 Then
                                                
                                                If (ExtremeTrim = 1 And BestXOList(A, B).Daughter = CP) Or BestXOList(A, B).MinorP = CP Or BestXOList(A, B).MajorP = CP Then
                                                    If ExtremeTrim = 1 Then
                                                        'ie remove this event if a an earlier triplet involved sequence x,y
                                                        TestC(SuperEventList(BestXOList(A, B).Eventnumber)) = 1
                                                        Y = YH
                                                        X = -1
                                                    Else
                                                        If BestXOList(A, B).MinorP = CP Then
                                                            ProbPar = BestXOList(A, B).MinorP
                                                        Else
                                                            ProbPar = BestXOList(A, B).MajorP
                                                        End If
                                                        
                                                        BPos2 = BestXOList(A, B).Beginning
                                                        EPos2 = BestXOList(A, B).Ending
                                                        If BPos2 < EPos2 Then
                                                            For z = BPos2 To EPos2
                                                                If FragReg(z) > 0 Then
                                                                    If AgeEvent(1, SuperEventList(BestXOList(A, B).Eventnumber)) < CurAge Then
                                                                        TestC(SuperEventList(BestXOList(A, B).Eventnumber)) = 1
                                                                        'A = PermNextNo
                                                                        'B = BCurrentXOver(PermNextNo)
                                                                        Y = YH
                                                                        X = -1
                                                                    
                                                                        Exit For
                                                                    End If
                                                                End If
                                                            Next z
                                                        Else
                                                            For z = 1 To EPos2
                                                                If FragReg(z) > 0 Then
                                                                    If AgeEvent(1, SuperEventList(BestXOList(A, B).Eventnumber)) < CurAge Then
                                                                        TestC(SuperEventList(BestXOList(A, B).Eventnumber)) = 1
                                                                        'A = PermNextNo
                                                                        'B = BCurrentXOver(PermNextNo)
                                                                        Y = YH
                                                                        X = -1
                                                                    
                                                                        Exit For
                                                                    End If
                                                                End If
                                                            Next z
                                                            If z = EPos + 1 Then
                                                                For z = BPos2 To Len(StrainSeq(0))
                                                                    If FragReg(z) > 0 Then
                                                                        If AgeEvent(1, SuperEventList(BestXOList(A, B).Eventnumber)) < CurAge Then
                                                                            TestC(SuperEventList(BestXOList(A, B).Eventnumber)) = 1
                                                                            'A = PermNextNo
                                                                            'B = BCurrentXOver(PermNextNo)
                                                                            Y = YH
                                                                            X = -1
                                                                        
                                                                            Exit For
                                                                        End If
                                                                    End If
                                                                Next z
                                                            End If
                                                        End If
                                                    End If
                                                End If
                                            End If
                                        End If
                                    'End If
                                End If
                            Next B
                        Next A
                    End If
                End If
            End If
        End If
    Next Y
Next X
If EarliestEvent <= SEventNumber Then
    Fail = EarliestEvent
End If
End Sub
Public Sub WriteTextNum(TP As PictureBox, S1 As String, S2)
Dim NumX As String, PowX As String, ExtraX As String, OS As String, OY As Long, OF As Double, NY As Long
Call NumToString(S2, 3, NumX, PowX, ExtraX)
OS = S1 + NumX + ExtraX

OY = TP.CurrentY
TP.CurrentX = 5
TP.Print OS
NY = TP.CurrentY
TP.CurrentY = OY
TP.CurrentX = 5 + TP.TextWidth(OS)
OF = TP.FontSize
TP.FontSize = 5
TP.Print PowX
TP.FontSize = OF
TP.CurrentY = NY
End Sub

Public Sub DoRDP(SPX As Byte, CPermNo As Long)
SAbortFlag = 0
MaxRepeatCycles = 2
Dim tDaught() As Long, OSX As String
DoQuick = 0
Dim tDon() As Byte
'ssss = GetTickCount
Dim AgeScore() As Double, EventScore() As Long, Y As Long, otMatch As Double, OKSeq() As Double, RFF As Byte, SRCompatF(2) As Double, SRCompatS(2) As Double, RList2() As Long, RNum2() As Long
Dim PhPrScore3() As Double, tPhPrScore3() As Double, SubPhPrScore3() As Double, SubScore3() As Double
Dim Consensus() As Double
Dim TDScores() As Double
Dim RetrimFlag As Byte
Dim PhPrScore2() As Double, tPhPrScore2() As Double, SubPhPrScore2() As Double, SubScore2() As Double
Dim SSDist(2) As Double
Dim oRNumX(2), UsedPar() As Long, SimScore() As Byte, SimScoreB() As Double, UNF() As Byte, INList() As Byte, OUList() As Byte, AAF As Byte, MIF As Byte, MaxBP(1) As Double, BPlots() As Double, oPMax As Double, oPMin As Double, RWinPP As Byte, NPh As Long, SetTot() As Long, TotD As Double, RCorrWarn() As Byte, DontRedo() As Byte
Dim BackRlist2() As Long, BackRNum2() As Long, BackRlist() As Long, BackRNum() As Long, JumpFlag As Byte, tEventAdd As Long, EPos As Long, BPos As Long, tRCorr() As Double, DMatS() As Double, LDst(1) As Double, AVSN(3) As Double, LowP As Double, RedoCycle As Byte, SLS As Long, UBx As Long
Dim InvS() As Byte, DoneProg() As Byte, TempDone() As Byte, DoneSeq() As Byte
Dim Outlyer(3) As Byte, SeqPair(2) As Byte, MinPair(1) As Byte
Dim DoneTarget As Byte, MissIDFlag As Byte, hMatch As Byte
Dim Relevant2() As Byte, PDist(2, 3) As Double, TMatch(1) As Double, TrpScore(3) As Double, PhPrScore() As Double, SubScore(3) As Double, SubPhPrScore(2) As Double
Dim MinDist(1) As Double, LDistXF(2) As Double, LDistXS(2) As Double, LDist(2) As Double, LDist2(2) As Double, LDist3(2) As Double, LDist4(2) As Double, LDistB(2) As Double, LDistB2(2) As Double, LDistB3(2) As Double, LDistB4(2) As Double, TtX(1) As Double, IntVal(1) As Double, MinDistZ(2) As Double
Dim BadDists() As Double, CoDists() As Double, ListCorr3() As Double, ListCorr() As Double, ListCorr2() As Double, tListCorr() As Double, LMat() As Double, RMat() As Double, tMat() As Double, RInv() As Double, RCorr() As Double, TDiffs() As Double, TValid() As Double, NumInGroup() As Double, MinPos() As Double, MaxPos() As Double, SubValidx() As Double, SubDiffsx() As Double, TotP() As Double, tVal() As Double, TPVal() As Double, PScores() As Double, RLScore() As Double
Dim LRC As Double, LPP As Double, HTS As Double, SNRD As Double
Dim MScore As Double, TPS As Double, CPVal As Double, BestPX As Double, TmF As Double, AvDst As Double

Dim CXO As Long, GoOnX As Long, ActualSeqSize() As Long, TraceSub() As Long, UnInvolved() As Long
Dim RCompatC(2) As Long, RCompatD(2) As Long, RCompatB(2) As Long, RCompat(2) As Long, BPosLR(3) As Long, RSize(5) As Long, CSeq(1) As Long, CompMat(2, 1) As Long, ISeqs(2) As Long, Trace(1) As Long
Dim NScoresX() As Long, WinnerPos() As Long, WinnerPosMa() As Long, WinnerPosMi() As Long, SLookUp() As Long, SLookUpNum() As Long, DoPairs() As Long, RListX() As Long, RNumX() As Long, RList() As Long, RNum() As Long, oRnum() As Long
Dim OuCheck() As Long, oBreaks() As Long, oRlist() As Long, tMatchX() As Long, FoundOne() As Long, SQ() As Long, OLSeq() As Long, OLSeqB() As Long, OLSeqE() As Long, GoodC() As Long, MinSeq() As Long, MaxSeq() As Long, GroupSeq() As Long, DoneOne() As Long, Breaks() As Long, InvListX() As Long, InvList() As Long, InPen() As Long, InPenX() As Long
Dim NumD As Long, X As Long, Mi As Long, DA As Long, Ma As Long, UB As Long
Dim OLSize As Long, TWinner As Long, BPos2 As Long, EPos2 As Long, tWinPP As Long, OS As Long, OE As Long, WinPP As Long, ActualE As Long, BXOSize As Long
Dim xNextno As Long, WinPPY As Long, EventAdd As Long, SP(5) As Long, EP(5) As Long, VSN As Long, TSN As Long, OldY As Long, XOSize As Long, WinRL As Long, oTotRecs As Long, SNextno As Long
Dim PhylCheck As Byte

Dim CAcList() As Integer, AcList() As Integer, EList() As Integer, tSeqNum() As Integer
Dim tBXOListMa() As XOverDefine, tBXOListMi() As XOverDefine, collectevents() As XOverDefine, CollectEventsMi() As XOverDefine, CollectEventsMa() As XOverDefine, PXOList() As XOverDefine, TempXOList() As XOverDefine, BestXOList() As XOverDefine
Dim tBCurrentXoverMa() As Integer, tBCurrentXoverMi() As Integer, DoneX() As Long, RCats() As Long, NRNum() As Long, NRList() As Long, NRNum2() As Long, NRList2() As Long, Relevant() As Long, PCurrentXOver() As Integer, TCurrentXOver() As Integer, BCurrentXOver() As Integer


oTotRecs = oRecombNo(100)

TotT = 0
For X = 0 To AddNum - 1
    TotT = TotT + MethodTime(X)
Next X

If TotT > 0 Then
    For X = 0 To AddNum - 1
        TimeFract(X) = MethodTime(X) / TotT
        X = X
    Next X
Else
    For X = 0 To AddNum - 1
        TimeFract(X) = 0.2
        
    Next X
    TotT = 1.2
End If
'reDim UsedPar(Len(StrainSeq(0)), NextNo)

If ShowPlotFlag = 2 And (CLine = "" Or CLine = " ") Then
        Call SetUpAxes(Nextno, XOverList(), CurrentXOver(), BPlots(), MaxBP(), oPMax, oPMin)
End If

'XX = doscans(0,2)
Form1.SSPanel1.Caption = "Looking For Unique Events"

LT = GetTickCount
TmF = 10 ^ -14
BXOSize = 10

Rnd (-3)
CHEvFlag = -1

RedoCycle = 0
ReDim SQ(2), RecombNo(AddNum), OLSeq(Len(StrainSeq(0))), OLSeqB(Len(StrainSeq(0))), OLSeqE(Len(StrainSeq(0))), Steps(4, 100)
ReDim ExtraHits(PermNextNo, 1), ExtraHitsMa(PermNextNo, 1), ExtraHitsMi(PermNextNo, 1)





If Nextno < 3 Then
    ForcePhylE = 0
End If
StepNo = 0
If DoQuick = 0 Then
    PhylCheck = 1
Else
    PhylCheck = 0
End If
'xxx = CurDir
'Open "Scores.csv" For Output As #1
'Close #1
oGCIndelFlag = GCIndelFlag
ReDim XDiffpos(Len(StrainSeq(0)) + 200), XPosDiff(Len(StrainSeq(0)) + 200)


    ReDim XoverSeqNum(Len(StrainSeq(0)), 2)
    ReDim ValidSpacer(Nextno), SpacerSeqs(Nextno)
    ReDim XOverSeqNumW(Len(StrainSeq(0)) + Int(XOverWindowX / 2) * 2, 2)



    'GCIndelFlag = 0
    If GCtripletflag = 1 Then
        ReDim SubSeq(Len(StrainSeq(0)), 6)
        ReDim FragMaxScore(Len(StrainSeq(0)), 5), MaxScorePos(Len(StrainSeq(0)), 5)
        ReDim PVals(Len(StrainSeq(0)), 5), FragSt(Len(StrainSeq(0)), 6), FragEn(Len(StrainSeq(0)), 6), FragScore(Len(StrainSeq(0)), 6)
    End If




    If MCTripletFlag = 0 Then
            HWindowWidth = CLng(MCWinSize / 2)
            ReDim Scores(Len(StrainSeq(0)), 2)  ' 0=s1,s2Matches etc
            ReDim Winscores(Len(StrainSeq(0)) + HWindowWidth * 2, 2) ' 0=s1,s2Matches etc
            ReDim ChiVals(Len(StrainSeq(0)), 2), ChiPvals(Len(StrainSeq(0)), 2), SmoothChi(Len(StrainSeq(0)), 2)
            If MCProportionFlag = 0 Then
                
                Call GetCriticalDiff(0)
                
                If MCWinSize <> HWindowWidth * 2 And MCProportionFlag = 0 Then
                    MCWinSize = HWindowWidth * 2
                End If
            End If
    End If


    HWindowWidth = CLng(CWinSize / 2)
    ReDim ScoresX(Len(StrainSeq(0)))  ' 0=s1,s2Matches etc
    ReDim WinScoresX(Len(StrainSeq(0)) + HWindowWidth * 2) ' 0=s1,s2Matches etc
    ReDim ChiValsX(Len(StrainSeq(0))), ChiPValsX(Len(StrainSeq(0))), SmoothChiX(Len(StrainSeq(0)))
    Call GetCriticalDiff(1)
        If CWinSize <> HWindowWidth * 2 And CProportionFlag = 0 Then
            CWinSize = HWindowWidth * 2
        End If
    


Dim WeightMod() As Long, Scratch() As Integer
ReDim Scratch(BSStepWin), WeightMod(BSBootReps, BSStepWin - 1)
Dummy = SEQBOOT2(BSRndNumSeed, BSBootReps, BSStepWin, Scratch(0), WeightMod(0, 0))


        Dim OnlySiScan As Byte
        If DoScans(0, 0) = 0 And DoScans(0, 1) = 0 And DoScans(0, 2) = 0 And DoScans(0, 3) = 0 And DoScans(0, 4) = 0 Then
            OnlySiScan = 1
        Else
            OnlySiScan = 0
        End If
        CurrentCorrect = 5
        Dim VRandTemplate() As Byte, HRandTemplate() As Long, TakenPos() As Byte
        
        'Dimension horizontal randomisation array if necessary
        If SSOutlyerFlag = 0 Or SSOutlyerFlag = 1 Then
            ReDim HRandTemplate(SSWinLen)
            ReDim TakenPos(SSWinLen)
        End If

        'Dimension vertical randomisation array
        ReDim VRandTemplate(Len(StrainSeq(0)), SSNumPerms)
        
        Dim DoGroupS() As Byte, DoGroupP() As Byte, DG1() As Byte, DG2() As Byte, VRandConv(15, 12) As Byte, Seq34Conv() As Byte
        ReDim DoGroupP(1, 3), DoGroupS(1, 3), DG1(15), DG2(14), Seq34Conv(5, 5)
        
        Call SetUpSiScan(Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP())
        LSeq = Len(StrainSeq(0))
        WinNum = CLng(Len(StrainSeq(0)) / SSStep + 1)
        Dim CorrectP As Double, TSub() As Long, oSeq As Long, PermSScores() As Long, PScoreHolder() As Long, SScoreHolder() As Long, PermPScores() As Long, SeqScore3() As Integer, MeanPScore() As Double, SDPScore() As Double
        ReDim SeqScore3(Len(StrainSeq(0))), MeanPScore(15), SDPScore(15)
        ReDim PermSScores(SSNumPerms, 15), PermPScores(SSNumPerms, 15), SScoreHolder(WinNum, 14), PScoreHolder(WinNum, 15)
        
        Dim SeqMap() As Byte, ZPScoreHolder() As Double, ZSScoreHolder() As Double
        ReDim SeqMap(Len(StrainSeq(0)))
        ReDim ZPScoreHolder(WinNum, 15)
        ReDim ZSScoreHolder(WinNum, 14)
        
        'ReDim PermPScores(SSNumPerms, 15)
        If SSOutlyerFlag = 2 Then
            Call GetOutie
            oSeq = Outie
        End If
        ReDim TSub(Nextno)
        For X = 0 To Nextno
            TSub(X) = X
        Next X
        
        
        If MCFlag = 0 Then
            CorrectP = LowestProb / MCCorrection
        Else
            CorrectP = LowestProb
        End If


ReDim Scores(Len(StrainSeq(0)), 2)

'xxxxzzzz ReDim EventSeq(3, Len(StrainSeq(0)), 1)

Eventnumber = 0

'Make permanent copy of xoverlist
XOSize = UBound(XOverList, 2)


ReDim TraceSub(Nextno), PermStraiName(Nextno)

For X = 0 To PermNextNo
    TraceSub(X) = X
    PermStraiName(X) = StraiName(TraceSub(X))
Next X
If Nextno > PermNextNo Then
    For X = PermNextNo + 1 To Nextno
        TraceSub(X) = S2TraceBack(X)
        PermStraiName(X) = StraiName(TraceSub(X))
    Next X
End If


ReDim SimSeqNum(Len(StrainSeq(0)), Nextno)


ReDim Preserve SeqNum(Len(StrainSeq(0)), Nextno), UsedPar(Len(StrainSeq(0)), Nextno)

TXOS = Nextno
'PermNextNo = NextNo
oNextno = Nextno
Dummy = CopySeqs(Len(StrainSeq(0)), Nextno, SeqNum(0, 0), SimSeqNum(0, 0))

CompMat(0, 0) = 1: CompMat(0, 1) = 2: CompMat(1, 0) = 0: CompMat(1, 1) = 2: CompMat(2, 0) = 0: CompMat(2, 1) = 1

STime = Abs(ST)

Call BuildFirstXOList(0, SPX, AgeScore(), EventScore(), MinSeqSize, JumpFlag, MissingData(), TraceSub(), Nextno, StepNo, Steps(), ExtraHits(), ExtraHitsMa(), ExtraHitsMi(), NOPINI(), Eventnumber, SEventNumber, BestXOList(), BCurrentXOver(), XOverList(), CurrentXOver(), BestXOListMi(), BCurrentXoverMi(), BestXOListMa(), BCurrentXoverMa(), Daught(), MinorPar(), MajorPar())

If JumpFlag = 0 Then
    ReDim BCurrentXoverMa(Nextno), BCurrentXoverMi(Nextno), BestXOListMa(Nextno, 10), BestXOListMi(Nextno, 10)
    ReDim BCurrentXOver(Nextno)
    ReDim BestXOList(Nextno, BXOSize), BestXOListMa(Nextno, BXOSize), BestXOListMi(Nextno, BXOSize)
End If

ReDim RepeatCycles(1000)
RestartX:
If SAbortFlag = 1 And AbortFlag = 0 Then
    AbortFlag = 1
    SAbortFlag = 0
End If
ReDim Relevant(Nextno), UnInvolved(Nextno)
ReDim PCurrentXOver(Nextno), TempXOList(Nextno, XOSize), PXOList(Nextno, XOSize)
ReDim DoneSeq(Nextno, UBound(PXOList, 2)), StepSEn(100)
ReDim NumRecsI(Nextno), Relevant2(2, Nextno)
MSize = 0
For X = 0 To Nextno
    For Y = 1 To CurrentXOver(X)
        If MSize < XOverList(X, Y).Daughter Then MSize = XOverList(X, Y).Daughter
        If MSize < XOverList(X, Y).MinorP Then MSize = XOverList(X, Y).MinorP
        If MSize < XOverList(X, Y).MajorP Then MSize = XOverList(X, Y).MajorP
    Next Y
Next X


Call CopyXOLists(XOSize, DoneSeq(), TempXOList(), PXOList(), PCurrentXOver(), XOverList(), CurrentXOver(), NumRecsI())

'Make permanent copy of seqnum


'Work out actual sequence sizes at the start
ReDim ActualSeqSize(Nextno)
For X = 0 To Nextno
    For Y = 1 To Len(StrainSeq(0))
        If SeqNum(Y, X) > 46 Then
            ActualSeqSize(X) = ActualSeqSize(X) + 1
        End If
    Next Y
Next X


'Open "numr.csv" For Output As #1
'    Close #1
'XX = SeqCol(3) ''3445930


Do 'loop until everything is completed from best to worst
    SNextno = Nextno
    DoneTarget = 0
    C = 0
    
    
    Do
       
        'find the best P-value
        'start Section 2
        LowP = LowestProb
        oRecombNo(100) = 0
        If X = X Then
            Do
                oRecombNo(100) = FindBestRecSignal(DoneTarget, Nextno, UBound(PXOList, 1), LowP, DoneSeq(0, 0), Trace(0), PCurrentXOver(0), PXOList(0, 0))
                
                If LowP = LowestProb Then Exit Do
                If PXOList(Trace(0), Trace(1)).MajorP = PXOList(Trace(0), Trace(1)).MinorP Or PXOList(Trace(0), Trace(1)).MajorP = PXOList(Trace(0), Trace(1)).Daughter Or PXOList(Trace(0), Trace(1)).Daughter = PXOList(Trace(0), Trace(1)).MinorP Then
                    LowP = LowestProb
                    PXOList(Trace(0), Trace(1)).Probability = 1
                    'XX = PXOList(Trace(0), Trace(1)).Eventnumber
                Else
                    Exit Do
                End If
            Loop
        Else
            Call FindBestRecSignalVB(oRecombNo(100), LowP, Trace(), PCurrentXOver(), DoneSeq(), PXOList())
        End If
        
        ET = Abs(GetTickCount)
        
        If Abs(ET - LT) > 500 Then
            LT = ET
            UpdateRecNums (SEventNumber)
            Form1.Label57(0).Caption = DoTimeII(Abs(ET - STime))
            If Abs(ET - elt) > 2000 Then
                elt = ET
                If oTotRecs > 0 Then
                    pbv = (1 - (oRecombNo(100) ^ 0.4 / oTotRecs ^ 0.4)) * 100
                    If pbv > Form1.ProgressBar1 Then
                        Form1.ProgressBar1 = pbv
                    End If
                End If
                                        
            End If
        End If
        
        If LowP = LowestProb And RedoListSize > 0 And RedoCycle = 0 Then
            If RedoCycle = 1 Then
                MaxXOListSize = MaxXOListSize * 2
                
            Else
                RedoCycle = 1
                
            End If
            ReDim CurrentXOver(Nextno)
            ReDim DonePVCO(AddNum - 1, Nextno), MaxXOP(AddNum - 1, Nextno)
            X = RedoListSize
            SLS = RedoListSize
            Call SignalCount(PXOList(), PCurrentXOver())
            Call UpdateRecNums(SEventNumber)
            Do While X >= 0
                GoOn = 1
                
                If GoOn = 1 Then
                    Seq1 = RedoList(1, X)
                    Seq2 = RedoList(2, X)
                    Seq3 = RedoList(3, X)
                    
                    If Seq1 <= Nextno And Seq2 <= Nextno And Seq3 <= Nextno Then
                       
                        If RedoList(0, X) = 0 Then
                            Call XOver(SeqNum(), Seq1, Seq2, Seq3, 0)
                        ElseIf RedoList(0, X) = 1 Then
                            Call GCXoverD(0)
                        ElseIf RedoList(0, X) = 2 Then
                             Call BSXoverS(BackUpNextno, Seq1, Seq2, Seq3, TraceSub(), MissingData())
                        ElseIf RedoList(0, X) = 3 Then
                            Call MCXoverF(0, 0, 0)
                            X = X
                        ElseIf RedoList(0, X) = 4 Then
                            tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                                                        
                            Call CXoverA(0, 0, 0)
                                                        
                            Seq1 = tseq2: Seq2 = tseq3: Seq3 = tseq1
                                                        
                            Call CXoverA(0, 0, 0)
                                                        
                            Seq1 = tseq3: Seq2 = tseq1: Seq3 = tseq2
                                                        
                            Call CXoverA(0, 0, 0)
                                                        
                            Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                        
                        ElseIf RedoList(0, X) = 5 Then
                            oSeq1 = Seq1: oSeq2 = Seq2: oSeq3 = Seq3
                            Call SSXoverC(0, WinNum, SeqMap(), ZPScoreHolder(), ZSScoreHolder(), CorrectP, oSeq, PermSScores(), PermPScores(), SScoreHolder(), PScoreHolder(), TraceSub(), SeqScore3(), MeanPScore(), SDPScore(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP())
                            Seq1 = oSeq1: Seq2 = oSeq2: Seq3 = oSeq3
                        ElseIf RedoList(0, X) = 6 Then
                        
                        End If
                    End If
                    If X >= 0 Then
                        RedoList(0, X) = -1
                    End If
                
                End If
                X = X - 1
                TT = GetTickCount
                If TT - LT > 500 Then
                    LT = TT
                    Form1.SSPanel1.Caption = Trim(Str(SLS - X)) & " of " & Trim(Str(SLS + 1)) & " triplets reexamined"
                End If
                
            Loop
            'clean up redolist
            Dummy = CleanRedoList(RedoListSize, RedoList(0, 0))
            
            Call CopyXOLists(XOSize, DoneSeq(), TempXOList(), PXOList(), PCurrentXOver(), XOverList(), CurrentXOver(), NumRecsI())
           
        
        ElseIf LowP = LowestProb Or AbortFlag = 1 Then
            If DoneTarget = 1 Or ForcePhylE = 1 Or AbortFlag = 1 Then
                Form1.ProgressBar1.Value = 100
                
                '2.797 old treedistances
                '2.766 new treedistances not in C
                '2.641 ints in RDP's Findseq and better nesting
                '2.328 ints in xohomol
                '2.265 not using pointers in xohomol
                '2.187 discalcx optimised quickdist
                '2.157 not using distance pointer in neighbor
                '2.141 better iterations in neighbour
                '1.863 using dopairs (but on laptop)
                '1.609 dopairs and optimised section 1.
                '1.578 - improvements in section 5
                '1.312 - improvements in section 9
                '1.250 - improvements in section 11/xover
                '1.142 (laptop) inprovements in section 6
                '0.941 (laptop) improvements in sections 16 - 19
                '0.871 (laptop) improvements in sections 10 and 12
                '0.801 (laptop)getting rid of section 9
                '0.781 (laptop) improvements in section 4
                '0.721 (laptop) improvements in section 6 & 7 using "relevant" array
                '0.500 (dsktop) even with using missing data
                '0.469 (dsktop) using checksplits & findmissing
                '0.590 (laptop)
                '0.581 -(laptop) maketrpgroups
                '0.551 - (laptop) makebposlr
                '0.422 (dsktop) -better dims in dordp
                '0.406 (dsktop) improvements in section 7
                '0.421
                '0.501 (laptop)
                '0.375 (dsktop)
                '0.359 using copyseq
                
                '7.515 no MC
                '7.020 (laptop no mc)
                '3.896
                '3.815 - makephprscore
                '3.605 - maketrpgroups
                '3.595 - makebposlr
                '3.956
                '3.825
                '3.535 - improvements in UPGMA
                '3.250 -improvements in dimentioning
                '3.047 - using makerlist
                '2.969 - using testrlist
                '3.365 (laptop)
                '2.735 (dsktop)
                '2.703 - makesplit
                '2.625 - stripunfound
                '2.610
                
                '1.522 (all evidence)
                '0.931 (forcephyle)
                '0.752 using highenough
                '0.721 calckmax
                '0.691 distributing hits accross xolist
                
                
                '2.584 (no correction)
                
                
                'maxchi results
                
                '0.681 laptop(correction)
                '0.541 lptop
                '0.531 laptop findhiseqs
                '0.511 using vquickdist
                
                '79.711 (no correction)
                '2.964 (laptop) wo quickdist
                '2.664 - using vquickdist
                '18.672 - maxxolistsize = 3 RGC no correct
                '17.750 - maxxolistsize = 6 RGC no correct
                '14.456 - maxxolistsize = 12 RGC no correct
                '9.578 - no maxxolistsize restriction
                
                '1.072
                '0.751
                '0.741
                '0.701
                Call Finalise(BCurrentXOver(), BestXOList())
                Form1.ProgressBar1.Value = 0
                Exit Sub
            Else
                DoneTarget = 1
            End If
        Else
            
             '1813,1871,5,9,10, 1.4^-2
            'Test movement in the tree
            
            EPos = PXOList(Trace(0), Trace(1)).Ending
            BPos = PXOList(Trace(0), Trace(1)).Beginning
            
            If BPos = 0 Then BPos = 1
            If EPos = 0 Then EPos = 1
            
            ISeqs(0) = PXOList(Trace(0), Trace(1)).Daughter
            ISeqs(1) = PXOList(Trace(0), Trace(1)).MinorP
            ISeqs(2) = PXOList(Trace(0), Trace(1)).MajorP
            
            
            LongFlag = DoneTarget
            'start Section 1***********************************************
            
            ReDim FAMat(Nextno, Nextno), SAMat(Nextno, Nextno)
            'SS = GetTickCount
            'For i = 0 To 10
            Call TestMoveInTree(1, BPos, EPos, SeqPair(), MinPair(), ISeqs(), SeqNum())
            'Next i
            'ee = GetTickCount
            'tt = ee - SS
            
            'X = X
            '16.944 total
            '9.834
            '8.913
            '5.608, 5.338, 5.408, 5.548
            '5.638, 5.828,5.408,5.358
            
            
            'XX = TraceSub(33)
            'end section 1*************************************************
            '12.749 5K perms
            '9.894 (9.814 without array dimentioning)
            '3.329 5K perms
            
            'if this doesn't move go back and find one that does
            
            If AbortFlag = 0 And (MinPair(0) <> MinPair(1) Or DoneTarget = 1) Then
                
                ' If SEventNumber >= 0 And X = X Then
                '    Open "numr.csv" For Append As #1
                '    Write #1,  oRecombNo(100)
                '    Close #1
              '
               ' End If
           
                For X = 0 To Nextno
                    FCMat(X, X) = 0
                    SCMat(X, X) = 0
                    
                Next X
                Form1.SSPanel1.Caption = "Finding daughter sequence"
                'is the diff in minpair translatable to a change in tree shape?
                MinDist(0) = 1000000
                MinDist(1) = 1000000
                
                Outlyer(0) = 2
                Outlyer(1) = 1
                Outlyer(2) = 0
                z = 0
                
               
                
                'Eventnumber = Eventnumber
                'Nextno = Nextno
                For X = 0 To 1
                    For Y = X + 1 To 2
                        If FAMat(ISeqs(X), ISeqs(Y)) < MinDist(0) Then
                            MinDist(0) = FAMat(ISeqs(X), ISeqs(Y))
                            MinPair(0) = z
                            SeqPair(0) = X
                            SeqPair(1) = Y
                            SeqPair(2) = Outlyer(z)
                        End If
                        If SAMat(ISeqs(X), ISeqs(Y)) < MinDist(1) Then
                            MinDist(1) = SAMat(ISeqs(X), ISeqs(Y))
                            MinPair(1) = z
                        End If
                        z = z + 1
                    Next Y
                Next X
                 
                RedoCycle = 0
                Eventnumber = Eventnumber + 1
                SEventNumber = SEventNumber + 1
                
                
                'Mark this specific example as the for this event.
                
                PXOList(Trace(0), Trace(1)).DHolder = (PXOList(Trace(0), Trace(1)).DHolder + 0.00000001) * -1
                 
                UB = UBound(StepSEn, 1)
                If SEventNumber > UB Then
                    ReDim Preserve StepSEn(UB + 100)
                End If
                
                
                Exit Do
                
            Else
                C = C + 1
                ET = Abs(GetTickCount)
                If Abs(ET - LT) > 500 Then
                    LT = ET
                    Form1.SSPanel1.Caption = Trim(Str(C)) & " of " & Trim(Str(oRecombNo(100))) & " rejected"
                    If Abs(ET - elt) > 2000 Then
                        elt = ET
                        If oTotRecs > 0 Then
                            pbv = (1 - (oRecombNo(100) ^ 0.4 / oTotRecs ^ 0.4)) * 100
                            If pbv > Form1.ProgressBar1 Then
                                Form1.ProgressBar1 = pbv
                            End If
                        End If
                                        
                    End If
                    UpdateRecNums (SEventNumber)
                    Form1.Label57(0).Caption = DoTimeII(Abs(ET - STime))
                End If
                
                DoneSeq(Trace(0), Trace(1)) = 1
            End If
        End If
    Loop
    
   
   'find smallest distance in famat
    
    
    'The best remaining event is found (favouring those that cause movement in the tree)
    'now test which is the recominant
    
    
    'begin section3******************
    
    'Call PhylProChecks(NextNo, SubPhPrScore(), SubScore(), PhPrScore(), SubValid(), SubDiffs(), PermValid(), PermDiffs(), FAMat(), SAMat(), FMat(), SMat(), Iseqs())
    
    Dim LD(1, 2, 1) As Double, TraceInvolvedBak() As Long, DoneThis() As Long, tPhPrScore(2) As Double
    ReDim TraceInvolvedBak(Nextno)
    ReDim DoneThis(1, Nextno)
   
   
    'make donethis
    
       
    Dummy = MakeDoneThis(Nextno, FMat(0, 0), FAMat(0, 0), SAMat(0, 0), LD(0, 0, 0), DoneThis(0, 0), ISeqs(0))
    
    For X = 0 To Nextno
        If DoneThis(0, X) = 0 Then
            If SubValid(ISeqs(0), X) > 0 And SubValid(ISeqs(1), X) > 0 And SubValid(ISeqs(2), X) > 0 Then
                If (SubDiffs(ISeqs(0), X) / SubValid(ISeqs(0), X)) > 0.6 Or SubDiffs(ISeqs(1), X) / SubValid(ISeqs(1), X) > 0.6 Or SubDiffs(ISeqs(2), X) / SubValid(ISeqs(2), X) > 0.6 Then
                    DoneThis(0, X) = 1: DoneThis(1, X) = 1
                End If
            End If
            If PermValid(ISeqs(0), X) > 0 And PermValid(ISeqs(1), X) > 0 And PermValid(ISeqs(2), X) > 0 Then
                If PermDiffs(ISeqs(0), X) / PermValid(ISeqs(0), X) > 0.6 Or PermDiffs(ISeqs(1), X) / PermValid(ISeqs(1), X) > 0.6 Or PermDiffs(ISeqs(2), X) / PermValid(ISeqs(2), X) > 0.6 Then
                    DoneThis(0, X) = 1: DoneThis(1, X) = 1
                End If
            End If
        End If
    Next X
   
    Call MakeSSDistB(ISeqs(), SSDist(), FCMat(), FAMat(), FMat(), SMat(), Nextno, DoneThis())

    'Phylpro correlation check - this uses all sites and not just variable sites
    'I must see if a variable sites version will work better
    ReDim PhPrScore(2)
        
    Dummy = MakePhPrScore(Nextno, TmF, TraceInvolvedBak(0), DoneThis(0, 0), ISeqs(0), PhPrScore(0), FMat(0, 0), SMat(0, 0), tPhPrScore(0), SubPhPrScore(0), SubScore(0))
    
    
    ReDim PhPrScore2(2), tPhPrScore2(2), SubScore2(3), SubPhPrScore2(2)
    
    'Call PhylProChecks(NextNo, SubPhPrScore2(), SubScore2(), PhPrScore2(), SubValid(), SubDiffs(), PermValid(), PermDiffs(), FAMat(), SAMat(), FAMat(), SAMat(), Iseqs())
    Dummy = MakePhPrScore(Nextno, TmF, TraceInvolvedBak(0), DoneThis(0, 0), ISeqs(0), PhPrScore2(0), FAMat(0, 0), SAMat(0, 0), tPhPrScore2(0), SubPhPrScore2(0), SubScore2(0))
    
    
    ReDim PhPrScore3(2), tPhPrScore3(2), SubScore3(3), SubPhPrScore3(2), DoneThis(1, Nextno)
   
    'make donethis
    'XX = Trace(1)
       
    Dummy = MakeDoneThis(Nextno, FMat(0, 0), FCMat(0, 0), SCMat(0, 0), LD(0, 0, 0), DoneThis(0, 0), ISeqs(0))
        
    
    For X = 0 To Nextno
        If DoneThis(0, X) = 0 Then
            If SubValid(ISeqs(0), X) > 0 And SubValid(ISeqs(1), X) > 0 And SubValid(ISeqs(2), X) > 0 Then
                If (SubDiffs(ISeqs(0), X) / SubValid(ISeqs(0), X)) > 0.6 Or SubDiffs(ISeqs(1), X) / SubValid(ISeqs(1), X) > 0.6 Or SubDiffs(ISeqs(2), X) / SubValid(ISeqs(2), X) > 0.6 Then
                    DoneThis(0, X) = 1
                    DoneThis(1, X) = 1
                End If
            End If
            If PermValid(ISeqs(0), X) > 0 And PermValid(ISeqs(1), X) > 0 And PermValid(ISeqs(2), X) > 0 Then
                If PermDiffs(ISeqs(0), X) / PermValid(ISeqs(0), X) > 0.6 Or PermDiffs(ISeqs(1), X) / PermValid(ISeqs(1), X) > 0.6 Or PermDiffs(ISeqs(2), X) / PermValid(ISeqs(2), X) > 0.6 Then
                    DoneThis(0, X) = 1
                    DoneThis(1, X) = 1
                End If
            End If
        End If
    Next X
    Dummy = MakePhPrScore(Nextno, TmF, TraceInvolvedBak(0), DoneThis(0, 0), ISeqs(0), PhPrScore3(0), FCMat(0, 0), SCMat(0, 0), tPhPrScore3(0), SubPhPrScore3(0), SubScore3(0))
    
    'end make phprscore
    '2.674
    '1.442 5k perms
    '0.040 5K perms
    
    '2.694 10K perms
    '0.090 10K perms
    'end section3******************************************
    
    'begin section 4**********************************************
    'More direct tree check for recombination
    'test for associations that have changed between the
    'famat and samat matrices for the sequence triplet in question
                
    'go down tree from closest to furthest away for each of the seqs
    '-Scores for all seqs on all branches of FAMat are averaged over the whole
    'branch
      
    
   For X = 0 To 2
                    
        ReDim DoneOne(Nextno), GroupSeq(Nextno), NumInGroup(Nextno)
        TrpScore(X) = 0
        Dummy = MakeTrpGroups(X, Nextno, NumInGroup(0), CompMat(0, 0), ISeqs(0), DoneOne(0), GroupSeq(0), MinDistZ(0), FAMat(0, 0))
        Dummy = MakeTrpScore(X, Nextno, FAMat(0, 0), SAMat(0, 0), TrpScore(0), NumInGroup(0), ISeqs(0), GroupSeq(0))
                    
        
    Next X
    
    
    'end section 4*************************************
    '0.110 5k perms
    '(1) find NO, PI and NI
    Call MakeINList(INList(), OUList(), MinPair())
        
    
    'begin section 5 ***************************************************
    'begin section 5.1 ***************************************************
    
    
    'Check for similar ddetected events with other sequences
    'do this first by looking accross each of the breakpoints and looking
    'for correlation between a potential recombinant sequenence and other sequences
    '(1) Each breakpoint is tested independantly using a simple distance scan
    'of sequences spanning the breakpoint
    '(2) Existing evidence of recombination is scanned for potential "co-recombinants"
    'identified in (1)
    
    '(1) independant testing of bps
    
    'Find VSN vriable sites on either side of two breakpoints
    
    TSN = 0: VSN = 60 '(XOverWindowX * 3) + mcwin
    
    Dummy = MakeBPosLR(VSN, TSN, Len(StrainSeq(0)), BPos, EPos, SeqNum(0, 0), ISeqs(0), BPosLR(0), AVSN(0))
    
    'end section 5.1 ****************************************************************
    '1.432 5K perms
    '0.040
    
    'begin section 5.2 ****************************************************************
    ReDim tRCorr(2, 2, 4, Nextno), DMatS(3, Nextno, Nextno), InvS(2, Nextno), RMat(2, 5, Nextno), tMat(Nextno, Nextno), RCorr(2, 2, Nextno), RInv(2, 2, Nextno), LMat(2, 3, 1, Nextno)
    'RInv is used to collect information on the polarity of the correlation - ie
    'whether, if you switch the parent groups in the triplet correlation scan across
    'the breakpoints you get a correlation or not
    OS = 0: OE = 0
    Dim SubDiffsY() As Double, SubValidY() As Double, DumMatA() As Double, DumMatB() As Double
    ReDim SubDiffsY(4, Nextno, Nextno), SubValidY(4, Nextno, Nextno), DumMatA(Nextno, Nextno), DumMatB(Nextno, Nextno)
    
    SP(0) = BPosLR(0): EP(0) = BPos - 1
    SP(1) = BPos: EP(1) = BPosLR(1)
    SP(2) = BPosLR(2): EP(2) = EPos
    SP(3) = EPos + 1: EP(3) = BPosLR(3)
    SP(4) = BPos: EP(4) = EPos
    For X = 0 To 3
        If SP(X) > Len(StrainSeq(0)) Then
            SP(X) = SP(X) - Len(StrainSeq(0))
        ElseIf SP(X) < 1 Then
            SP(X) = SP(X) + Len(StrainSeq(0))
        End If
            
        If EP(X) > Len(StrainSeq(0)) Then
            EP(X) = EP(X) - Len(StrainSeq(0))
        ElseIf EP(X) < 1 Then
            EP(X) = EP(X) + Len(StrainSeq(0))
        End If
    Next X
    
    ReDim RCorrWarn(2)
   'XX = RCorrWarn(2)
    Call MakeProperRCorr(INList(), MissingData(), RCorrWarn(), RInv(), tRCorr(), RCorr(), CompMat(), Len(StrainSeq(0)), ISeqs(), SP(), EP(), SeqNum())
    If RCorrWarn(0) = 1 And RCorrWarn(1) = 1 Then RCorrWarn(2) = 0
   
    For X = 0 To 3
        
        If OS <> SP(X) Or OE <> EP(X) Then
            UB = UBound(PermValid, 1)
            ReDim SubValidx(Nextno, Nextno), SubDiffsx(Nextno, Nextno), DumMatB(Nextno, Nextno)
            'change this to only look at distances between all and iseqs(x)? done
            'Dummy = QuickDist2(Len(StrainSeq(0)), NextNo, UB, SP, EP, SubValidx(0, 0), SubDiffsx(0, 0), SeqNum(0, 0), ISeqs(0))
            Dummy = QuickDist(Len(StrainSeq(0)), Nextno, UB, SP(X), EP(X), DumMatA(0, 0), DumMatB(0, 0), SubValidx(0, 0), SubDiffsx(0, 0), PermValid(0, 0), PermDiffs(0, 0), SeqNum(0, 0))

        End If
       
        For Y = 0 To Nextno
            For z = Y + 1 To Nextno
                If SubValidx(Y, z) > 10 Then
                    DMatS(X, z, Y) = CLng(DumMatB(z, Y) * 10000) / 10000
                    DMatS(X, Y, z) = DMatS(X, z, Y)
                    SubValidY(X, Y, z) = SubValidx(Y, z)
                    SubValidY(X, z, Y) = SubValidx(Y, z)
                    SubDiffsY(X, Y, z) = SubDiffsx(Y, z)
                    SubDiffsY(X, z, Y) = SubDiffsx(Y, z)
                Else
                    DMatS(X, z, Y) = 3
                    DMatS(X, Y, z) = 3
                    SubValidY(X, Y, z) = 0
                    SubValidY(X, z, Y) = 0
                    SubDiffsY(X, Y, z) = 0
                    SubDiffsY(X, z, Y) = 0
                End If
            Next z
        Next Y
        
        OS = SP(X): OE = EP(X)
        
    Next X
    
    'check for possible rcorr problems
    
    'If X = X Then
        For X = 0 To 3
        
            TotD = DMatS(X, ISeqs(0), ISeqs(1)) + DMatS(X, ISeqs(0), ISeqs(2)) + DMatS(X, ISeqs(2), ISeqs(1))
            If TotD > 0 Then
                TotD = 2 / TotD
                If (1 - DMatS(X, ISeqs(0), ISeqs(1)) * TotD) < 0.4 And (1 - DMatS(X, ISeqs(0), ISeqs(2)) * TotD) < 0.4 And (1 - DMatS(X, ISeqs(2), ISeqs(1)) * TotD) < 0.4 Then
                    If X < 2 Then
                        RCorrWarn(0) = 1
                    Else
                        RCorrWarn(1) = 1
                    End If
                End If
            Else
                TotD = 0
                If X < 2 Then
                    RCorrWarn(0) = 1
                Else
                    RCorrWarn(1) = 1
                End If
                'make a warning
            End If
        Next X
        If RCorrWarn(0) = 1 And RCorrWarn(1) = 1 Then
            RCorrWarn(2) = 0
        ElseIf RCorrWarn(0) = 1 Or RCorrWarn(1) = 1 Then
            RCorrWarn(2) = 1
        End If
    'End If
    
   
    '****Test removing parental domination of R scores
    'PDist(0, 0) = (DMatS(0, ISeqs(1), ISeqs(2)) + DMatS(1, ISeqs(1), ISeqs(2))) / 2
    'PDist(0, 1) = PDist(0, 0)
    'PDist(0, 2) = (DMatS(2, ISeqs(1), ISeqs(2)) + DMatS(3, ISeqs(1), ISeqs(2))) / 2
    'PDist(0, 3) = PDist(0, 2)
   '
   ' PDist(1, 0) = (DMatS(0, ISeqs(0), ISeqs(2)) + DMatS(1, ISeqs(0), ISeqs(2))) / 2
   ' PDist(1, 1) = PDist(1, 0)
   ' PDist(1, 2) = (DMatS(2, ISeqs(0), ISeqs(2)) + DMatS(3, ISeqs(0), ISeqs(2))) / 2
   ' PDist(1, 3) = PDist(1, 2)
   '
   ' PDist(2, 0) = (DMatS(0, ISeqs(1), ISeqs(0)) + DMatS(1, ISeqs(1), ISeqs(0))) / 2
    ''PDist(2, 1) = PDist(2, 0)
    'PDist(2, 2) = (DMatS(2, ISeqs(1), ISeqs(0)) + DMatS(3, ISeqs(1), ISeqs(0))) / 2
    'PDist(2, 3) = PDist(2, 2)
    '**************************
    '1.610 up to here with 5K perms
    'end section 5 ****************************************************************
    '36.752 5K perms
    'end section 5.2 ****************************************************************
    '34.245 5K perms
    '32.719 (desktop)
    '32.500 (redims moved outside the loops)
    '6.531 '1K perms
    '5.812 - removed  superfluous loops
    '5.734 - ""
    '5.203 - don't bother multiplying everything by 1000
    '4.969 - replace 10^-14 with constants
    '4.532 - only calculating distances when necessary
    '4.343 quickdist2 instead of quickdist
    '3.297 using makermat
    '3.105 5K perms
    '2.143 5K perms
    '2.174 (laptop)
    
    'Find out how many steps needed if each of Iseqs were recombinant and
    'put into rcompat.
    'Do this by
    '(1) collecting all R values greater that 0.83 (r2 = 0.7)
    '(on either side of bp)
    '(2) finding most distant common ancestor of the group.
    '(3) Find offspring of the MRCA
    '(4) subtract potential recombinants from the offspring group
    '(5) go through the recombinant group and count the number of distance
    'categories between each recombinant and non-recombinants: rcompat = the
    'highest number of categories.
    
    'begin section 6 ***************************************************************
    'begin section 6.1 ***************************************************************
    
   
    'ss = GetTickCount
    'For i = 0 To 5000
    ReDim GoodC(Nextno, 1), InPenX(2), InPen(2), NRList(2, Nextno), NRNum(2), NRList2(2, Nextno), NRNum2(2), oRnum(2), RLScore(2, Nextno), RListX(2, Nextno, 9), RNumX(2, 9), InvListX(2, Nextno), RList(2, Nextno), RNum(2), InvList(2, Nextno), tVal(2, Nextno, 2), TPVal(2, Nextno, 2), PScores(2, Nextno, 2), TotP(2, Nextno, 1)
    
    'Find valid comparisons - anything with 10 or more overlapping nts is valid
    Dummy = MakeGoodC(Nextno, Len(StrainSeq(0)), GoodC(0, 0), BPosLR(0), SeqNum(0, 0))
    
    
    

   
   ' XX = GoodC(11, 0): XX = GoodC(11, 1)
   ' X = X
    '10.615 5K perms
    '0.187 5K perms
    
  
   
    
    
    'TotPTarget = 1
    'I must port this to C++ when I'm done debugging
    
    
    

    'Take another look at ssdist
    
    Dim OUIndexA(2) As Byte
    If SSDist(INList(0)) > SSDist(INList(1)) And SSDist(INList(0)) > SSDist(INList(2)) Then
        OUIndexA(INList(0)) = 1: OUIndexA(INList(1)) = 0: OUIndexA(INList(2)) = 0
    ElseIf SSDist(INList(0)) < SSDist(INList(1)) And SSDist(INList(0)) < SSDist(INList(2)) Then
        OUIndexA(INList(0)) = 0: OUIndexA(INList(1)) = 1: OUIndexA(INList(2)) = 1
    Else
        OUIndexA(INList(0)) = 0: OUIndexA(INList(1)) = 0: OUIndexA(INList(2)) = 0
    End If
    Dim AcceptableCoR() As Byte
   
    ReDim AcceptableCoR(2, Nextno)
    For X = 0 To Nextno
       
        If FAMat(ISeqs(INList(0)), X) < FAMat(ISeqs(INList(2)), ISeqs(INList(0))) Or SAMat(ISeqs(INList(1)), X) < SAMat(ISeqs(INList(1)), ISeqs(INList(0))) Then
            AcceptableCoR(INList(0), X) = 1 '0,1,2,3,4,9,10,11,12,13,14,16,18,19
            AcceptableCoR(INList(1), X) = 1
        End If
        
        If FAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(2)), ISeqs(INList(0))) Or SAMat(ISeqs(INList(2)), X) < SAMat(ISeqs(INList(2)), ISeqs(INList(0))) Then
            AcceptableCoR(INList(2), X) = 1 '0,1,2,3,12,13,15,18,19
        End If
    Next X
    
    Call MakeDontRedo(Nextno, ISeqs(), RList(), RNum(), FAMat(), SAMat(), INList(), RCorr(), DontRedo())
    
    
    For X = 0 To 2
    '5,1,5
        'RNumX(X) = RNum(X)
        For Y = 0 To RNum(X)
        
            For z = 0 To 9
                RListX(X, Y, z) = RList(X, Y)
            Next z
        Next Y
    Next X
    

    Dummy = MakeRList(Nextno, GoodC(0, 0), ISeqs(0), RListX(0, 0, 0), InvListX(0, 0), RNumX(0, 0), RList(0, 0), InvList(0, 0), RNum(0), RInv(0, 0, 0), RCorr(0, 0, 0), PScores(0, 0, 0), TPVal(0, 0, 0), tVal(0, 0, 0), TotP(0, 0, 0), RLScore(0, 0), DontRedo(0, 0), AcceptableCoR(0, 0), RCorrWarn(0))
    
    
    oRNumX(0) = RNum(0)
    oRNumX(1) = RNum(1)
    oRNumX(2) = RNum(2)
  ' XX = TPVal(2, 12, 0)
    ReDim OKSeq(2, 17, Nextno)
    Call AddOK1(PermNextNo, DontRedo(), RCorr(), RList(), RNum(), OKSeq(), TPVal(), RInv())
    
    
    'if there are contradictions go with the no inversion option.
   
    For Y = 0 To 2
        For X = 0 To RNum(Y)
        
            If InvList(Y, X) = 1 Then
                For z = 0 To 2
                    If RCorr(Y, z, RList(Y, X)) > 0.83 And RInv(Y, z, RList(Y, X)) = 0 Then
                        InvList(Y, X) = 0 '4,0,3:1,0,3:
                        For A = 0 To 2
                            If RInv(Y, A, RList(Y, X)) > 0 Then
                                RCorr(Y, A, RList(Y, X)) = 0 And RInv(Y, A, RList(Y, X)) = 0
                            End If
                        Next A
                    End If
                Next z
            End If
        Next X
    Next Y


    
    
    
    If MinPair(0) <> MinPair(1) Then
        
        
        
        'Make ListCorr
        'this is where I compare rcorrs and rinvs for evidence that each of the iseqs is recombinant
        '(1) find NO, PI and NI
        '(2) find expected lists for NO, PI and NI: NO=0, PI=1, NI=2
        '(3) for each rcor (left, right and middle or 0,1,2) work out which list fits the observed
        'correlations best
        
        

        Call MakeOUCheck(Nextno, ISeqs(), MinPair(), INList(), FAMat(), SAMat(), OuCheck())
        'Call MakeOUCheck(NextNo, ISeqs(), MinPair(), INList(), FCMat(), SCMat(), OuCheckB(), OuCheck2B(), OuCheck3B())
        
        
        
         '(2) find expected lists for NO, PI and NI: NO=0, PI=1, NI=2
         ReDim EList(2, 2, Nextno)
         For X = 0 To 2
             
             For Y = 0 To 2
                 For z = 0 To Nextno
                     EList(X, Y, z) = -1
                 Next z
             Next Y
         Next X
         'DO NO - ie we do the situation where the outlyer in the recombinant region is recombinant
         'we are assuming that iseqs(inlist(0)) is the recombinant and trying to
         'figure out what rlists would be obtained if this were the case
         
     
         
         For X = 0 To Nextno
             If FAMat(ISeqs(INList(0)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(1))) And SAMat(ISeqs(INList(0)), X) < SAMat(ISeqs(INList(0)), ISeqs(INList(1))) Then
                 EList(0, INList(0), X) = 0
             End If
             
             If FAMat(ISeqs(INList(0)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(2))) And FAMat(ISeqs(INList(0)), X) > 0 And SAMat(ISeqs(INList(1)), X) < SAMat(ISeqs(INList(0)), ISeqs(INList(1))) Then
                 EList(0, INList(1), X) = 0
             End If
             
             If FAMat(ISeqs(INList(0)), X) > FAMat(ISeqs(INList(0)), ISeqs(INList(2))) And SAMat(ISeqs(INList(0)), X) > SAMat(ISeqs(INList(2)), X) Then
                 EList(0, INList(1), X) = 2
             End If
             
             If FAMat(ISeqs(INList(0)), X) > FAMat(ISeqs(INList(0)), ISeqs(INList(1))) And SAMat(ISeqs(INList(0)), X) > SAMat(ISeqs(INList(1)), X) Then
                 EList(0, INList(2), X) = 0
             End If
             If FAMat(ISeqs(INList(0)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(1))) And FAMat(ISeqs(INList(0)), X) > 0 And SAMat(ISeqs(INList(0)), X) > SAMat(ISeqs(INList(1)), X) Then
                 EList(0, INList(2), X) = 2
             End If
             
         Next X
         
         'DO PI - ie we do the situation where the inlyer in both the recombinant region and
         'the background is recombinant
         'we are assuming that iseqs(inlist(1)) is the recombinant and trying to
         'figure out what rlists would be obtained if this were the case
         For X = 0 To Nextno
             If FAMat(ISeqs(INList(1)), X) > 0 And FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(2))) And SAMat(ISeqs(INList(1)), X) > SAMat(ISeqs(INList(0)), X) Then
                 EList(1, INList(0), X) = 0
             End If
             If FAMat(ISeqs(INList(2)), X) > 0 And FAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(2))) And SAMat(ISeqs(INList(2)), X) > SAMat(ISeqs(INList(1)), ISeqs(INList(2))) Then
                 EList(1, INList(0), X) = 1
             End If
             
             If FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(0))) And SAMat(ISeqs(INList(1)), X) < SAMat(ISeqs(INList(1)), ISeqs(INList(0))) Then
                 EList(1, INList(1), X) = 0
             End If
             
             If FAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(2))) And SAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(2))) Then
                 EList(1, INList(2), X) = 0
             End If
             If FAMat(ISeqs(INList(1)), X) > FAMat(ISeqs(INList(1)), ISeqs(INList(0))) And FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(2))) Then
                 EList(1, INList(2), X) = 1
             End If
             If FAMat(ISeqs(INList(1)), X) > 0 And FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(0))) And SAMat(ISeqs(INList(1)), X) < SAMat(ISeqs(INList(1)), ISeqs(INList(0))) Then
                 EList(1, INList(2), X) = 4
             End If
         Next X
       
         'DO NI - ie we do the situation where the outlyer in the background region is recombinant
         'we are assuming that iseqs(inlist(2)) is the recombinant and trying to
         'figure out what rlists would be obtained if this were the case
         For X = 0 To Nextno
             If FAMat(ISeqs(INList(1)), X) > 0 And FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(2))) Then
                 EList(2, INList(0), X) = 0
             End If
             If FAMat(ISeqs(INList(2)), X) > 0 And FAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(2))) And SAMat(ISeqs(INList(2)), X) < SAMat(ISeqs(INList(0)), X) Then
                 EList(2, INList(0), X) = 1
             End If
             
             If FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(0))) And SAMat(ISeqs(INList(1)), X) < SAMat(ISeqs(INList(1)), ISeqs(INList(0))) Then
                 EList(2, INList(1), X) = 0
             End If
             If FAMat(ISeqs(INList(2)), X) > 0 And FAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(2))) And SAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(0)), X) Then
                 EList(2, INList(1), X) = 4
             End If
             If FAMat(ISeqs(INList(1)), X) > FAMat(ISeqs(INList(1)), ISeqs(INList(0))) And FAMat(ISeqs(INList(1)), X) < FAMat(ISeqs(INList(1)), ISeqs(INList(2))) Then
                 EList(2, INList(1), X) = 2
             End If
             
             If FAMat(ISeqs(INList(2)), X) < FAMat(ISeqs(INList(0)), ISeqs(INList(2))) And SAMat(ISeqs(INList(2)), X) < SAMat(ISeqs(INList(0)), ISeqs(INList(2))) Then
                 EList(2, INList(2), X) = 0
             End If
             
         Next X
         
         '(3) for each rcor (left, right and middle or 0,1,2) work out which list fits the observed
         'correlations best
         'first build AcList
         ReDim AcList(2, 2, Nextno)
         'note that the 2nd dimention here refers to the beginning, end and middle
         'recombinant region
         
         For X = 0 To 2
             For z = 0 To 2
                 For Y = 0 To Nextno
                     If RInv(X, z, Y) = 3 Then RInv(X, z, Y) = 2
                 Next Y
             Next z
         Next X
         
         
        
         ReDim ListCorr3(2), ListCorr2(2), ListCorr(2), tListCorr(2, 2)
                  
         Call MakeListCorr(ISeqs(), INList(), EList(), tRCorr(), RNum(), RList(), RInv(), AcList(), ListCorr3(), ListCorr2(), ListCorr(), tListCorr(), RCorrWarn())
         'Call MakeListCorr(ISeqs(), INList(), EList(), tRCorr(), RNum(), RList(), RInv(), AcList(), ListCorr3(), ListCorr2(), ListCorr(), tListCorr())
    Else
        ReDim ListCorr3(2), ListCorr2(2), ListCorr(2), tListCorr(2, 2)
        ReDim OuCheck(2)

    
    End If
    
    
    For X = 0 To 2
        For Y = 0 To RNum(X)
            If InvList(X, Y) = 1 Then
                InvS(X, RList(X, Y)) = 1
                
            End If
        Next Y
    Next X
    
    
    
    
    'end of section 6.1******************************************************
    '11.486 5K perms
    '4.734 5K perms
    '4.281 - just taking out soem of the crap in ttest
    '1.797 - ttest using ttestprob
    '0.407 - using makerlist
    'check to see if events suggested in RList are there or not
    '2.047!!!!!
    '1.985
    
    '1.873 - 5kperms with 350 events
    
    'Set up array for region overlap
    
    
      
    
    Dummy = MakeOLSeq(Len(StrainSeq(0)), SP(0), EP(1), RSize(2), OLSeqB(0))
    Dummy = MakeOLSeq(Len(StrainSeq(0)), SP(2), EP(3), RSize(4), OLSeqE(0))
    Dummy = MakeOLSeq(Len(StrainSeq(0)), BPos, EPos, RSize(0), OLSeq(0))
    
    Dummy = MakeRelevant(Nextno, Relevant(0), RNum(0), RList(0, 0))
    
    Dim BMatch() As Double, BPMatch() As Long
    ReDim BMatch(2, Nextno), BPMatch(2, 1, Nextno)
    'BMatch(0) = 0: BMatch(1) = 0: BMatch(2) = 0
    '1.412
    ReDim SQ(3), UNF(2, Nextno)
    
    If IndividualA = -1 Then    'only do this if a proper balanced search for recombination
        
            
        If RedoListSize > 0 Then  'Need to see if any events in redolist need to be repeated
            ReDim CurrentXOver(Nextno), DonePVCO(AddNum - 1, Nextno), MaxXOP(AddNum - 1, Nextno)
            Call SignalCount(PXOList(), PCurrentXOver())
            Call UpdateRecNums(SEventNumber)
            For X = 0 To UBound(Relevant2, 2)
                Relevant2(0, X) = 0: Relevant2(1, X) = 0: Relevant2(2, X) = 0
            Next X
            For X = 0 To 2
                For Y = 0 To RNum(X)
                    Relevant2(X, RList(X, Y)) = 1
                Next Y
            Next X
            X = RedoListSize
            SLS = RedoListSize
            Do While X >= 0
                
                Dummy = FindNextRedo(X, Relevant2(0, 0), RedoList(0, 0), ISeqs(0), Relevant(0))
                If Dummy > 0 Then
                    X = Dummy
                    Seq1 = RedoList(1, X): Seq2 = RedoList(2, X): Seq3 = RedoList(3, X)
                    If PermValid(Seq1, Seq2) > MinSeqSize And PermValid(Seq1, Seq3) > MinSeqSize And PermValid(Seq2, Seq3) > MinSeqSize Then
                        If SubValid(Seq1, Seq2) > 20 And SubValid(Seq1, Seq3) > 20 And SubValid(Seq2, Seq3) > 20 Then
                            If Seq1 <= Nextno And Seq2 <= Nextno And Seq3 <= Nextno Then
                                If RedoList(0, X) = 0 Then
                                    Call XOver(SeqNum(), Seq1, Seq2, Seq3, 0)
                                ElseIf RedoList(0, X) = 1 Then
                                    Call GCXoverD(0)
                                ElseIf RedoList(0, X) = 2 Then
                                    Call BSXoverS(BackUpNextno, Seq1, Seq2, Seq3, TraceSub(), MissingData())
                                ElseIf RedoList(0, X) = 3 Then
                                    Call MCXoverF(0, 0, 0)
                                    X = X
                                ElseIf RedoList(0, X) = 4 Then
                                    tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                                                                
                                    Call CXoverA(0, 0, 0)
                                                                
                                    Seq1 = tseq2: Seq2 = tseq3: Seq3 = tseq1
                                                                
                                    Call CXoverA(0, 0, 0)
                                                                
                                    Seq1 = tseq3: Seq2 = tseq1: Seq3 = tseq2
                                                                
                                    Call CXoverA(0, 0, 0)
                                                                
                                    Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                                ElseIf RedoList(0, X) = 5 Then
                                    oSeq1 = Seq1: oSeq2 = Seq2: oSeq3 = Seq3
                                    Call SSXoverC(0, WinNum, SeqMap(), ZPScoreHolder(), ZSScoreHolder(), CorrectP, oSeq, PermSScores(), PermPScores(), SScoreHolder(), PScoreHolder(), TraceSub(), SeqScore3(), MeanPScore(), SDPScore(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP())
                                    Seq1 = oSeq1: Seq2 = oSeq2: Seq3 = oSeq3
                                ElseIf RedoList(0, X) = 6 Then
                                
                                End If
                            End If
                        End If
                    End If
                    If X >= 0 Then
                        RedoList(0, X) = -1
                    End If
                Else
                    Exit Do
                End If
                X = X - 1
                TT = GetTickCount
                If TT - LT > 500 Then
                
                    LT = TT
                    Form1.SSPanel1.Caption = Trim(Str(SLS - X)) & " of " & Trim(Str(SLS + 1)) & " triplets reexamined"
                End If
                
            Loop
            'clean up redolist
            Dummy = CleanRedoList(RedoListSize, RedoList(0, 0))
            
            
            'Add new events to pxolist
            Call CopyXOLists(XOSize, DoneSeq(), TempXOList(), PXOList(), PCurrentXOver(), XOverList(), CurrentXOver(), NumRecsI())
            Call SignalCount(PXOList(), PCurrentXOver())
            Call UpdateRecNums(SEventNumber)
        End If
        
       
    
        Call FindActualEventsVB(RLScore(), UNF(), InvList(), Nextno, RSize(), BPMatch(), BMatch(), OKSeq(), FoundOne(), SP(), EP(), RCorr(), OLSeq(), OLSeqB(), OLSeqE(), CSeq(), RNum(), RList(), InvS(), TMatch(), PXOList(), PCurrentXOver(), SQ(), tDon(), ISeqs(), CompMat())
    
        
    Else
        For X = 0 To 2
            For Y = 0 To RNum(X)
                UNF(X, Y) = 1
            Next Y
        Next X
    End If
    
    ReDim RList2(2, Nextno), RNum2(2)
    
    Call FindSets(OKSeq(), SetTot(), Nextno, BPos, EPos, ISeqs(), RList2(), RNum2(), PXOList(), PCurrentXOver())
    X = X
    '4.531
    '3.635
    '4.246
    
    '24.132
    '20.150 using relevant
    'strip duplicates and inversions from rlist
    
    Dummy = StripDupInv(Nextno, RCorr(0, 0, 0), RLScore(0, 0), InPen(0), RNum(0), RList(0, 0), InvList(0, 0))
    
    For X = 0 To 2
        For Y = 0 To RNum(X)
            OKSeq(X, 4, RList(X, Y)) = 1
        Next Y
    Next X
    
    '******************need to add in a similarity check for high corrolaters*******************
    '******************only accept if corrolaters more similer over BP than  *******************
    '******************the rest of the sequences in other groups are to one  *******************
    '******************another                                               *******************
    
 '   For X = 0 To 2
 '       Y = 0
 '       Do While Y <= RNumX(X)
 '           For Z = 0 To RNum(X)
 '               If RListX(X, Y) = RList(X, Z) Then
 '                   Exit For
 '               End If
 '           Next Z
 '           If Z > RNum(X) Then
 '               If Y < RNumX(X) Then
 '                   RListX(X, Y) = RListX(X, RNumX(X))
 '               End If
 '               RNumX(X) = RNumX(X) - 1
 '           Else
 '               Y = Y + 1
 '           End If
 '       Loop
 '       If InPenX(X) = 1 And InPen(X) = 0 Then InPenX(X) = 0
 '   Next X
    Dim RCompatNF(2, 9), RCompatNS(2, 9)
     
    For X = 0 To 2
        For Y = 0 To 9
            RCompatNF(X, Y) = 0: RCompatNS(X, Y) = 0
        Next Y
    Next X
    
    Dim RCompatXF(2) As Long, RCompatXS(2) As Long, RCompatBXF(2) As Long, RCompatBXS(2) As Long, RCompat2(2) As Long, RCompatB2(2) As Long, RCompat3(2) As Long, RCompatB3(2) As Long, RCompat4(2) As Long, RCompatB4(2) As Long
    Dim RCompatS(2) As Long, RCompatBS(2) As Long, RCompatS2(2) As Long, RCompatBS2(2) As Long, RCompatS3(2) As Long, RCompatBS3(2) As Long, RCompatS4(2) As Long, RCompatBS4(2) As Long
        
   For WinPP = 0 To 2
        
        RCompat(WinPP) = 0: RCompatB(WinPP) = 0: RCompat2(WinPP) = 0: RCompatB2(WinPP) = 0
        RCompat3(WinPP) = 0: RCompatB3(WinPP) = 0: RCompat4(WinPP) = 0: RCompatB4(WinPP) = 0
        
        RCompatS(WinPP) = 0: RCompatBS(WinPP) = 0: RCompatS2(WinPP) = 0: RCompatBS2(WinPP) = 0
        RCompatS3(WinPP) = 0: RCompatBS3(WinPP) = 0: RCompatS4(WinPP) = 0: RCompatBS4(WinPP) = 0
    Next WinPP
    
    Dummy = MakeLDist(Nextno, LDist(0), FAMat(0, 0), RNum(0), RList(0, 0))
    Dummy = MakeLDist(Nextno, LDist3(0), SAMat(0, 0), RNum(0), RList(0, 0))
    
    'Dummy = MakeLDist(NextNo, LDistXF(0), FAMat(0, 0), RNumX(0), RListX(0, 0))
    'Dummy = MakeLDist(NextNo, LDistXS(0), SAMat(0, 0), RNumX(0), RListX(0, 0))
    
    For WinPP = 0 To 2
        ReDim DoneX(Nextno), RCats(Nextno * 3)
        Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompat(0), RCompatB(0), InPen(0), RCats(0), RNum(0), NRNum(0), GoodC(0, 0), DoneX(0), RList(0, 0), NRList(0, 0), FAMat(0, 0), LDist(0))
        'Call MakeRCompatVB(ISeqs(), Compmat(), WinPP, NextNo, RCompat(), RCompatB(), InPen(), RCats(), RNum(), NRNum(), GoodC(), DoneX(), RList(), NRList(), FAMat(), LDist())
        
        'Using SMat
        ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
        Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompatS(0), RCompatBS(0), InPen(0), RCats(0), RNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList(0, 0), NRList2(0, 0), SAMat(0, 0), LDist3(0))
        
        
        'Call MakeRCompatVB(ISeqs(), Compmat(), WinPP, NextNo, RCompatS(), RCompatBS(), InPen(), RCats(), RNum(), NRNum2(), GoodC(), DoneX(), RList(), NRList2(), SAMat(), LDist3())
        
        ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
       ' Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, NextNo, RCompatXF(0), RCompatBXF(0), InPenX(0), RCats(0), RNumX(0), NRNum2(0), GoodC(0, 0), DoneX(0), RListX(0, 0), NRList2(0, 0), FAMat(0, 0), LDistXF(0))
        ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
       ' Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, NextNo, RCompatXS(0), RCompatBXS(0), InPenX(0), RCats(0), RNumX(0), NRNum2(0), GoodC(0, 0), DoneX(0), RListX(0, 0), NRList2(0, 0), SAMat(0, 0), LDistXS(0))
        
      
            
    Next WinPP
    If RCompat(0) <> RCompat(1) Or RCompat(0) <> RCompat(2) Then
        
    Else
        If Nextno > 2 Then
            Dummy = MakeLDist(Nextno, LDist2(0), FCMat(0, 0), RNum(0), RList(0, 0))
            For WinPP = 0 To 2
                ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
                Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompat2(0), RCompatB2(0), InPen(0), RCats(0), RNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList(0, 0), NRList2(0, 0), FCMat(0, 0), LDist2(0))
            Next WinPP
        End If
        
        If RCompat2(0) <> RCompat2(1) Or RCompat2(0) <> RCompat2(2) Then
        Else
            Dummy = MakeLDist(Nextno, LDistB(0), FAMat(0, 0), RNum2(0), RList2(0, 0))
            For WinPP = 0 To 2
                'setsRcompat with non-BS tree
                ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
                Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompat3(0), RCompatB3(0), InPen(0), RCats(0), RNum2(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList2(0, 0), NRList2(0, 0), FAMat(0, 0), LDistB(0))
            Next WinPP
            
            If RCompat3(0) <> RCompat3(1) Or RCompat3(0) <> RCompat3(2) Then
            Else
                If Nextno > 2 Then
                    Dummy = MakeLDist(Nextno, LDistB2(0), FCMat(0, 0), RNum2(0), RList2(0, 0))
                    For WinPP = 0 To 2
                        ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
                        Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompat4(0), RCompatB4(0), InPen(0), RCats(0), RNum2(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList2(0, 0), NRList2(0, 0), FCMat(0, 0), LDistB2(0))
                    Next WinPP
                End If
            End If
        End If
    End If
    
    If RCompatS(0) <> RCompatS(1) Or RCompatS(0) <> RCompatS(2) Then
    
    Else
        If Nextno > 2 Then
            Dummy = MakeLDist(Nextno, LDist4(0), SCMat(0, 0), RNum(0), RList(0, 0))
            For WinPP = 0 To 2
                'normal rcompats with BS tree
                ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
                Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompatS2(0), RCompatBS2(0), InPen(0), RCats(0), RNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList(0, 0), NRList2(0, 0), SCMat(0, 0), LDist4(0))
            Next WinPP
        End If
        If RCompatS2(0) <> RCompatS2(1) Or RCompatS2(0) <> RCompatS2(2) Then
        Else
            Dummy = MakeLDist(Nextno, LDistB3(0), SAMat(0, 0), RNum2(0), RList2(0, 0))
            For WinPP = 0 To 2
                'setsrcompats with non-BS tree
                ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
                Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompatS3(0), RCompatBS3(0), InPen(0), RCats(0), RNum2(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList2(0, 0), NRList2(0, 0), SAMat(0, 0), LDistB3(0))
            Next WinPP
            If RCompatS3(0) <> RCompatS3(1) Or RCompatS3(0) <> RCompatS3(2) Then
            Else
                If Nextno > 2 Then
                    Dummy = MakeLDist(Nextno, LDistB4(0), SCMat(0, 0), RNum2(0), RList2(0, 0))
                    'setsrcompats with BS tree
                    For WinPP = 0 To 2
                        ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
                        Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompatS4(0), RCompatBS4(0), InPen(0), RCats(0), RNum2(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList2(0, 0), NRList2(0, 0), SCMat(0, 0), LDistB4(0))
                    Next WinPP
                End If
            End If
        End If
    End If
   ' For X = 0 To 2
   '         For Y = 0 To RNum(X)
   '         'XX = InvList(X, Y)
  '          'XX = ISeqs(X)
  '          'XX = StraiName(ISeqs(X))
  '          'XX = StraiName(RList(X, Y))
   '             XX = RList(X, Y) '19,14,11
   '         Next Y
   '     Next X
    
    'Clean up rlistx
    Dim FEntry() As Byte
    'For X = 0 To 2
    '    For Y = 0 To 9
    '        XX = RListX(2, 0, 0) 'RNumX(1, Y)
    '        XX = RNum(X)
    '    Next Y
    'Next X
    
    For X = 0 To 2
    'XX = ISeqs(0)
        ReDim FEntry(Nextno)
        For Y = 0 To RNum(X)
            FEntry(RList(X, Y)) = 1
        Next Y
        For Y = 1 To 9 'ie leave 0 alone
            z = 0
            Do While z <= RNumX(X, Y)
                If FEntry(RListX(X, z, Y)) = 0 Then
                    If z < RNumX(X, Y) Then
                        RListX(X, z, Y) = RListX(X, RNumX(X, Y), Y)
                    End If
                    RNumX(X, Y) = RNumX(X, Y) - 1
                Else
                    z = z + 1
                End If
            Loop
            
        Next Y
    Next X
    
    Dim txLDist(2) As Double, txRNum(2) As Long, txRList() As Long
    For X = 0 To 9
        For WinPP = 0 To 2
            RCompatXF(WinPP) = 0: RCompatXS(WinPP) = 0: RCompatBXF(WinPP) = 0: RCompatBXS(WinPP) = 0
        Next WinPP
        
        txLDist(0) = 0: txLDist(1) = 0: txLDist(2) = 0
        ReDim txRList(2, Nextno)
        
        For Y = 0 To 2
            txRNum(Y) = RNumX(Y, X)
            For z = 0 To RNumX(Y, X)
                txRList(Y, z) = RListX(Y, z, X)
            Next z
        Next Y
        
        Dummy = MakeLDist(Nextno, txLDist(0), FAMat(0, 0), txRNum(0), txRList(0, 0))
        For WinPP = 0 To 2
            ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno), InPenX(2)
            Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompatXF(0), RCompatBXF(0), InPenX(0), RCats(0), txRNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), txRList(0, 0), NRList2(0, 0), FAMat(0, 0), txLDist(0))
        Next WinPP
        
        RCompatNF(0, X) = RCompatXF(0): RCompatNF(1, X) = RCompatXF(1): RCompatNF(2, X) = RCompatXF(2)
        
        txLDist(0) = 0: txLDist(1) = 0: txLDist(2) = 0
        ReDim txRList(2, Nextno)
        
        For Y = 0 To 2
            txRNum(Y) = RNumX(Y, X)
            For z = 0 To RNumX(Y, X)
                txRList(Y, z) = RListX(Y, z, X)
            Next z
        Next Y
        
        Dummy = MakeLDist(Nextno, txLDist(0), SAMat(0, 0), txRNum(0), txRList(0, 0))
        For WinPP = 0 To 2
            ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno), InPenX(2)
            Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), WinPP, Nextno, RCompatXS(0), RCompatBXS(0), InPenX(0), RCats(0), txRNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), txRList(0, 0), NRList2(0, 0), SAMat(0, 0), txLDist(0))
        Next WinPP
        
        RCompatNS(0, X) = RCompatXS(0): RCompatNS(1, X) = RCompatXS(1): RCompatNS(2, X) = RCompatXS(2)
    Next X
    
    For Y = 0 To 2
        SRCompatF(Y) = 0: SRCompatS(Y) = 0
        For X = 0 To 9
            SRCompatF(Y) = SRCompatF(Y) + RCompatNF(Y, X)
            SRCompatS(Y) = SRCompatS(Y) + RCompatNS(Y, X)
        Next X
        SRCompatF(Y) = SRCompatF(Y) / 10 + InPen(Y)
        SRCompatS(Y) = SRCompatS(Y) / 10 + InPen(Y)
        
    Next Y
    
    For WinPP = 0 To 2
        RCompatXF(WinPP) = 0: RCompatXS(WinPP) = 0: RCompatBXF(WinPP) = 0: RCompatBXS(WinPP) = 0
    Next WinPP
    
    
    
    
    ReDim SimScore(2)
    Call SimpleDist(Nextno, SimScore(), SimScoreB(), RList(), RNum(), INList(), ISeqs(), FMat(), SMat())
    ReDim CoDists(2)
    
    'work out baddists
    Call GetBadDists(Nextno, RCorr(), FAMat(), CompMat(), DMatS(), ISeqs(), BadDists(), UNF(), RList(), RNum())
   
    
    'draw trees with regions centred on the breakpoints
    
    If PhylCheck = 1 Then
            Dim BPMat() As Double
            Call MakeBPMatX(Nextno, BPMat(), FMat(), DMatS(), GoodC())
    End If
    
    
    
    ReDim Consensus(2, 3), TDScores(20, 2)
        
    'Call MakeConsensus(TDScores(), Consensus(), INList(), CompMat(), SubScore2(), SubPhPrScore2(), OuCheck2(), OuCheck3(), ListCorr3(), SubScore(), PhPrScore2(), ListCorr(), ListCorr2(), TrpScore(), PhPrScore(), SubPhPrScore(), OuCheck(), BadDists(), CoDists(), RCompat(), tListCorr())
    
    
    RetrimFlag = 0
    If RCompat(0) > 0 And RCompat(1) > 0 And RCompat(2) > 0 And RCompatS(0) > 0 And RCompatS(1) > 0 And RCompatS(2) > 0 Then
        RetrimFlag = 1
    End If
    For X = 0 To 2
        RCompatC(X) = 0: RCompatD(X) = 0
    Next X
     
    'XX = RNum(0)
    '    XX = RNum(1)
    '    XX = RNum(2)
    '    RNum(0) = 0
    '    RNum(1) = 0
    '    RNum(2) = 0
    '    RList(0, 0) = ISeqs(0)
    '    RList(1, 0) = ISeqs(1)
    '    RList(2, 0) = ISeqs(2)
    '    X = X
  'Start trimming sequences out of Rlist
  
    Call CheckPattern(OKSeq(), ISeqs(), RList(), SP(), EP(), CompMat(), RNum(), SeqNum(), RInv(), InvList(), RCorrWarn())
     
           
    MissIDFlag = 0
    RFF = 0
    If RetrimFlag = 1 Then
        RWinPP = 4
        RFF = 1
        Call FinalTrim(WinNum, SeqMap(), ZPScoreHolder(), ZSScoreHolder(), LowP, 0, WeightMod(), OKSeq(), BMatch(), BPMatch(), MinPair(), INList(), RInv(), BPos, EPos, RWinPP, MissIDFlag, PhylCheck, BackUpNextno, Nextno, CorrectP, oSeq, PermSScores(), PermPScores(), SScoreHolder(), PScoreHolder(), SeqScore3(), MeanPScore(), SDPScore(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP(), XOverList(), CurrentXOver(), PCurrentXOver(), PXOList(), MissingData(), TraceSub(), Trace(), DMatS(), CompMat(), BPMat(), RCorrWarn(), UNF(), RCorr(), ISeqs(), FCMat(), SCMat(), FAMat(), SAMat(), NRNum(), NRList(), RNum(), RList())
        
        For X = 0 To 2
            RCompatB(X) = 0
            ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
            Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), X, Nextno, RCompatC(0), RCompatB(0), InPen(0), RCats(0), RNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList(0, 0), NRList2(0, 0), FAMat(0, 0), LDist(0))
            RCompatB(X) = 0
            ReDim DoneX(Nextno), RCats(Nextno * 3), NRNum2(2), NRList2(2, Nextno)
            Dummy = MakeRCompat(ISeqs(0), CompMat(0, 0), X, Nextno, RCompatD(0), RCompatB(0), InPen(0), RCats(0), RNum(0), NRNum2(0), GoodC(0, 0), DoneX(0), RList(0, 0), NRList2(0, 0), SAMat(0, 0), LDist3(0))
        'do RCompatC
        Next X
    End If
    
    'Dim ConsensusB() As Double
    'ReDim ConsensusB(2, 2)
    
    'Call MakeConsensusB(RCompatC(), RCompatD(), RCompat2(), RCompat3(), RCompat4(), RCompatS(), RCompatS2(), RCompatS3(), RCompatS4(), SSDist(), OUIndexA(), TDScores(), ConsensusB(), INList(), CompMat(), SubScore2(), SubPhPrScore2(), ListCorr3(), SubScore(), PhPrScore2(), ListCorr(), ListCorr2(), TrpScore(), PhPrScore(), SubPhPrScore(), OuCheck(), BadDists(), CoDists(), RCompat(), tListCorr())
    Call MakeConsensusC(RCompatXF(), RCompatXS(), SimScore(), SimScoreB(), PhPrScore3(), RCompatC(), RCompatD(), RCompat2(), RCompat3(), RCompat4(), RCompatS(), RCompatS2(), RCompatS3(), RCompatS4(), SSDist(), OUIndexA(), TDScores(), Consensus(), INList(), CompMat(), SubScore2(), SubPhPrScore2(), ListCorr3(), SubScore(), PhPrScore2(), ListCorr(), ListCorr2(), TrpScore(), PhPrScore(), SubPhPrScore(), OuCheck(), BadDists(), CoDists(), RCompat(), tListCorr())



    
    If X = 12345 Then
        
        If SEventNumber = 1 Then
            Open "Scores.csv" For Output As #1
            Close #1
        End If
        
        Open "scores.csv" For Append As #1
        Write #1, SEventNumber
        Dim ccon(2) As Byte
        
        Dim NNumb() As Long
        ReDim NNumb(2)
        For A = 0 To 2
            For C = 1 To Len(StraiName(ISeqs(A)))
                NNumb(A) = NNumb(A) + Asc(Mid(StraiName(ISeqs(A)), C, 1)) * 2 ^ (Len(StraiName(ISeqs(A))) - C)
            Next C
        Next
        
        'B,A,03:1252,1202,1271
        'b,03,a:1252,1271,1202
        If NNumb(0) < NNumb(1) And NNumb(0) < NNumb(2) Then
            ccon(0) = 0
            If NNumb(1) < NNumb(2) Then
                ccon(1) = 1
                ccon(2) = 2
            Else
                ccon(2) = 1
                ccon(1) = 2
            End If
        ElseIf NNumb(1) < NNumb(2) Then
            ccon(0) = 1
            If NNumb(0) < NNumb(2) Then
                ccon(1) = 0
                ccon(2) = 2
            Else
                ccon(2) = 0
                ccon(1) = 2
            End If
        
        Else
            ccon(0) = 2
            If NNumb(0) < NNumb(1) Then
                ccon(1) = 0
                ccon(2) = 1
            Else
                ccon(2) = 0
                ccon(1) = 1
            End If
            X = X
        End If
        For A = 0 To 2
        
            C = ccon(A)
            
            If PermNextNo > 2 Then
                Write #1, StraiName(TraceSub(ISeqs(C))), oRNumX(C), RNum(C) + 1, ListCorr(C), NPh, SimScoreB(C), SimScore(C), PhPrScore(C), PhPrScore2(C), PhPrScore3(C), SubScore(C), SSDist(C), OUIndexA(C), SubPhPrScore(C), SubScore2(C), SubPhPrScore2(C), SRCompatF(C), SRCompatS(C), RCompat(C), RCompat2(C), RCompat3(C), RCompat4(C), RCompatS(C), RCompatS2(C), RCompatS3(C), RCompatS4(C), RCompatXF(C), RCompatXS(C), RCompatC(C), RCompatD(C), TrpScore(C), BadDists(C), OUList(C), ListCorr2(C), ListCorr3(C), Consensus(C, 0), Consensus(C, 1), Consensus(C, 2), OuCheck(C), SetTot(0, C), SetTot(1, C)    ', PhPrScore3(C), SubScore3(C), SubPhPrScore3(C),   TrpScore2(C), RCompat2(C) ', SubScore2(C), SubPhPrScore2(C), ListCorr4(C), tListCorr2(C, 0), tListCorr2(C, 1), tListCorr2(C, 2), ListCorr5(C), ListCorr6(C)
            Else
                Write #1, StraiName(TraceSub(ISeqs(C))), RNum(C) + 1, ListCorr(C), NPh, PhPrScore(C), PhPrScore2(C), SubScore(C), SubPhPrScore(C), SubScore2(C), SubPhPrScore2(C), RCompat(C), RCompat2(C), RCompat3(C), RCompat4(C), RCompatS(C), RCompatS2(C), RCompatS3(C), RCompatS4(C), TrpScore(C), BadDists(C), OUList(C), ListCorr2(C), ListCorr3(C), Consensus(C, 0), Consensus(C, 1), Consensus(C, 2), OuCheck(C), SetTot(0, C), SetTot(1, C)   ', PhPrScore3(C), SubScore3(C), SubPhPrScore3(C), OuCheckB(C),  TrpScore2(C), RCompat2(C) ', SubScore2(C), SubPhPrScore2(C), ListCorr4(C), tListCorr2(C, 0), tListCorr2(C, 1), tListCorr2(C, 2), ListCorr5(C), ListCorr6(C)
            End If
        Next A
        Write #1, ""
        BPos = BPos
        EPos = EPos
        Close #1
    End If
    
    
    
    'Now guestimate which are winpp (ie the whole list if necessary)
    
    'First check to see if one of the sequences is very poorly aligned with the others
    
    For X = 0 To 2
        If Consensus(X, 2) >= Consensus(CompMat(X, 0), 2) And Consensus(X, 2) >= Consensus(CompMat(X, 1), 2) Then
            
            WinPP = X
        End If
    Next X
    
    If (Consensus(WinPP, 2) / (Consensus(0, 2) + Consensus(1, 2) + Consensus(2, 2))) < 0.6 Then
        MissIDFlag = MissIDFlag + 10
    End If
    Call QTestAlign(WinPP, LDst(), ISeqs(), PermDiffs(), PermValid(), SubDiffs(), SubValid())
    
    AAF = 0
    'MissIDFlag = 0
    
    If LDst(1) > 0.6 And CPermNo = 0 And X = 12345 Then
        MissIDFlag = MissIDFlag + 3
        For X = 0 To 2
            RNum(X) = 0
            RList(X, 0) = ISeqs(X)
        Next X
        Call DelSeq(BPos, EPos, ISeqs(WinPP))
        AAF = 1
    Else
        
        If CPermNo = 0 And Realignflag = 1 And LDst(1) > 0.2 Then
            Call RAlignAndRecheck(PXOList(Trace(0), Trace(1)).ProgramFlag, SeqNum(), ISeqs(), BPos, EPos, AAF)
            If AAF <> 1 Then
                zzztX = zzztX + 1
            Else
                zzzt = zzzt + 1
            End If
        End If
        If AAF = 1 Then
            
            MissIDFlag = MissIDFlag + 3
            For X = 0 To 2
                RNum(X) = 0
                RList(X, 0) = ISeqs(X)
            Next X
            
            Call DelSeq(BPos, EPos, ISeqs(WinPP))
        Else
            
            
             If LDst(1) > 0.5 Or AAF = 2 Then
                 MissIDFlag = MissIDFlag + 2
                 
             Else
                 MissIDFlag = MissIDFlag + 0
             End If
             
             'find lowest rcompat
            
             
             
        End If
    End If
    Dim TBreak(2) As Double
    For X = 0 To 2
        TBreak(X) = PhPrScore(X) - TrpScore(X) + RCompat(X) + RCompatS(X)
    Next X
    If Consensus(WinPP, 2) = Consensus(CompMat(WinPP, 0), 2) And Consensus(WinPP, 2) = Consensus(CompMat(WinPP, 1), 2) Then
        If TBreak(0) < TBreak(1) And TBreak(0) < TBreak(2) Then
            WinPP = 0
        ElseIf TBreak(1) < TBreak(0) And TBreak(1) < TBreak(2) Then
            WinPP = 1
        Else
            WinPP = 2
        End If
    ElseIf Consensus(WinPP, 2) = Consensus(CompMat(WinPP, 0), 2) Then
        If TBreak(WinPP) > TBreak(CompMat(WinPP, 0)) Then
            WinPP = CompMat(WinPP, 0)
        
        End If
    ElseIf Consensus(WinPP, 2) = Consensus(CompMat(WinPP, 1), 2) Then
        If TBreak(WinPP) > TBreak(CompMat(WinPP, 1)) Then
            WinPP = CompMat(WinPP, 1)
        
        End If
    End If
    
    
        
        tWinPP = WinPP
        Dim tListX() As Byte
        ReDim NRList(2, Nextno)
        For WinPP = 0 To 2
        
            
            ReDim tListX(Nextno)
            
            For X = 0 To RNum(WinPP)
                tListX(RList(WinPP, X)) = 1
            Next X
            NRNum(WinPP) = 0
            For X = 0 To Nextno
                If tListX(X) = 0 Then
                    NRList(WinPP, NRNum(WinPP)) = X
                    NRNum(WinPP) = NRNum(WinPP) + 1
                End If
            Next X
            NRNum(WinPP) = NRNum(WinPP) - 1
        Next WinPP
        WinPP = tWinPP
    
  
    'edit RList to contain only those sequences that are not separeted by
    'non-recombinants from the original hit - ie iseqs(winpp)
    'This is probably where the identity check should be included....
    
    'Add parents to nrlist
   ' For X = 0 To NRNum(WinPP)
   '     XX = NRList(WinPP, X)
   ' Next X
    
    'For X = 0 To 2
    '    If X <> WinPP Then
    '        nrlnum = NRNum(WinPP) + 1
    '
    '    End If
    'Next X
   
'************* Got to expand NRList to include all sequences not in RList**************
    
    If AAF <> 1 Then
        'If SEventNumber = 19 Then
        '    XX = TraceSub(50)
        'End If
        RWinPP = WinPP
        Call FinalTrim(WinNum, SeqMap(), ZPScoreHolder(), ZSScoreHolder(), LowP, RFF, WeightMod(), OKSeq(), BMatch(), BPMatch(), MinPair(), INList(), RInv(), BPos, EPos, RWinPP, MissIDFlag, PhylCheck, BackUpNextno, Nextno, CorrectP, oSeq, PermSScores(), PermPScores(), SScoreHolder(), PScoreHolder(), SeqScore3(), MeanPScore(), SDPScore(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP(), XOverList(), CurrentXOver(), PCurrentXOver(), PXOList(), MissingData(), TraceSub(), Trace(), DMatS(), CompMat(), BPMat(), RCorrWarn(), UNF(), RCorr(), ISeqs(), FCMat(), SCMat(), FAMat(), SAMat(), NRNum(), NRList(), RNum(), RList())
        'Call FinalTrim2(RCorr(), MinPair(), FCMat(), SCMat(), INList(), ISeqs(), RNum(), RList())
        WinPP = RWinPP
    End If
    
    'Keep a note of which sequences have been used as parents and where
     
     'For X = 0 To RNum(0)
     '   XX = RList(0, X) '19,31,32,33,34,37,38,40,46,47,48,50,52,54,55,56,57,58,59,60,61,62,63,65
     'Next X
    'Call CheckUsedPar(WinPP, SEventNumber, BPos, EPos, Compmat(), UsedPar(), MissingData())
    
    
  ' XX = TraceSub(62)
    
    
    'begin section 7*****************************************************
    'Find events most closely resembling the RList events in PXOList
    '(1) scan for best P-value for events with identical parents and method
    '(2) scan for "worse" events for other methods but with highest P-value
    'For X = 0 To 2
    '    BackRNum2(X) = RNum(X)
    '    For Y = 0 To RNum(X)
    '        BackRlist2(X, Y) = RList(X, Y)
    '    Next Y
    'Next X
    Dummy = MakeRelevant(Nextno, Relevant(0), RNum(0), RList(0, 0))
    
    '6.360
    ReDim collectevents(Nextno, AddNum), CollectEventsMa(Nextno, AddNum), CollectEventsMi(Nextno, AddNum)
    
    If AAF <> 1 Then
        SS = GetTickCount
        Call MakeCollecteventsX(SMat(), ISeqs(), LowP, WinPP, EPos, BPos, OLSeq(), RSize(), CompMat(), Relevant(), PCurrentXOver(), collectevents(), PXOList(), RNum(), RList())
    
        Call MakeCollecteventsX(SMat(), ISeqs(), LowP, CompMat(WinPP, 0), EPos, BPos, OLSeq(), RSize(), CompMat(), Relevant(), PCurrentXOver(), CollectEventsMa(), PXOList(), RNum(), RList())
        
        Call MakeCollecteventsX(SMat(), ISeqs(), LowP, CompMat(WinPP, 1), EPos, BPos, OLSeq(), RSize(), CompMat(), Relevant(), PCurrentXOver(), CollectEventsMi(), PXOList(), RNum(), RList())
        EE = GetTickCount
    TT = EE - SS
     X = X
    Else
        collectevents(0, PXOList(Trace(0), Trace(1)).ProgramFlag) = PXOList(Trace(0), Trace(1))
        CollectEventsMa(0, PXOList(Trace(0), Trace(1)).ProgramFlag) = PXOList(Trace(0), Trace(1))
        CollectEventsMi(0, PXOList(Trace(0), Trace(1)).ProgramFlag) = PXOList(Trace(0), Trace(1))
    End If
   
     
    'end section 7********************************************
    '4.727 5K perms
    '3.835 using olsize
    '2.864 with relevant
    '2.375 (dsktop)
    '1.750 with maketmatch2
    '1.7641 - taking out unnecessary crap
    '1.062 - with tetsrlist
   ' Next i
    
    'ee = GetTickCount
    'tt = ee - ss
    'X = X
    '22.762 10 times more events
    '21.892
    '14.331 with relevant
    '9.109 (dsktop) with testrlist
    
    'begin section 8*****************************************************
    
    
    ReDim tSeqNum(Len(StrainSeq(0)), RNum(WinPP)), WinnerPos(RNum(WinPP), AddNum), WinnerPosMa(RNum(CompMat(WinPP, 0)), AddNum), WinnerPosMi(RNum(CompMat(WinPP, 1)), AddNum)
    UB1 = UBound(ExtraHits, 1)
    UB2 = UBound(ExtraHitsMa, 1)
    UB3 = UBound(ExtraHitsMi, 1)
    ReDim Preserve ExtraHits(UB1, SEventNumber), ExtraHitsMa(UB2, SEventNumber), ExtraHitsMi(UB3, SEventNumber)
    
    If ShowPlotFlag = 2 And CLine = "" Then
        
        Dim DistX As Double
        If MissIDFlag = 13 Or MissIDFlag = 3 Then
            Mi = 1
        Else
            Mi = 0
            Call UpdatePlotsF(BPlots(), BPos, EPos, MaxBP(), PXOList(Trace(0), Trace(1)).SBPFlag)
        End If
        
        DistX = 10
        For X = 0 To 2
            For Y = X + 1 To 2
                If FMat(ISeqs(X), ISeqs(Y)) < SMat(ISeqs(X), ISeqs(Y)) Then
                    If DistX > FMat(ISeqs(X), ISeqs(Y)) Then
                        DistX = FMat(ISeqs(X), ISeqs(Y))
                    End If
                Else
                    If DistX > SMat(ISeqs(X), ISeqs(Y)) Then
                        DistX = SMat(ISeqs(X), ISeqs(Y))
                    End If
                End If
            Next Y
        Next X
        Call UpdatePlotE(BPos, EPos, PXOList(Trace(0), Trace(1)).Probability, oPMax, oPMin, DistX, Mi, MaxBP(), BPlots())
    End If
    
    Call MakeBestXOList(EventAdd, BXOSize, WinPP, RCorr(), ExtraHits(), SuperEventList(), TraceSub(), WinnerPos(), RList(), RNum(), BCurrentXOver(), BestXOList(), collectevents())
    
    
    
    
    
    ReDim tBXOListMa(Nextno, AddNum - 1), tBXOListMi(Nextno, AddNum - 1), tBCurrentXoverMa(Nextno), tBCurrentXoverMi(Nextno)
    'XX = SEventNumber
    Dim tEHitsMa() As Byte, tEHitsMi() As Byte
    ReDim tEHitsMa(Nextno), tEHitsMi(Nextno)
    Call MakeNextBestXOLists(RCorr(), tEHitsMa(), tEventAdd, CompMat(WinPP, 0), SuperEventList(), TraceSub(), WinnerPosMa(), RList(), RNum(), tBCurrentXoverMa(), tBXOListMa(), CollectEventsMa())
    If tEventAdd > EventAdd Then EventAdd = tEventAdd
    Call MakeNextBestXOLists(RCorr(), tEHitsMi(), tEventAdd, CompMat(WinPP, 1), SuperEventList(), TraceSub(), WinnerPosMi(), RList(), RNum(), tBCurrentXoverMi(), tBXOListMi(), CollectEventsMi())
    If tEventAdd > EventAdd Then EventAdd = tEventAdd
     
    'end section 8********************************************
    '9.994 5K perms
    '9.794 - incorrect because of SEventnumber going too high
    '0.503 = actual time
    
    'save info on the exact sequences used
    
    '10.194 10 times more events
    
    'begin section 9*******************************************
'xxxxzzzz    ReDim Preserve EventSeq(3, Len(StrainSeq(0)), Eventnumber + EventAdd)
    
    
    Dim SwapFlag As Byte
    ReDim Preserve NOPINI(2, SEventNumber)
    Dim DMiMa(2) As Long
    If INList(0) <> INList(1) And INList(0) <> INList(2) And INList(1) <> INList(2) Then
        If INList(0) = WinPP Then 'NO recombinant
            NOPINI(0, SEventNumber) = 0: NOPINI(1, SEventNumber) = 1: NOPINI(2, SEventNumber) = 2
            DMiMa(INList(0)) = 0: DMiMa(INList(1)) = 1: DMiMa(INList(2)) = 2
            
            If INList(1) = CompMat(WinPP, 0) Then
                SwapFlag = 0
            Else
                SwapFlag = 1
            End If
        ElseIf INList(1) = WinPP Then 'PI recombinant
            
            NOPINI(0, SEventNumber) = 1: NOPINI(1, SEventNumber) = 0: NOPINI(2, SEventNumber) = 2
            DMiMa(INList(0)) = 1: DMiMa(INList(1)) = 0: DMiMa(INList(2)) = 2
            If INList(0) = CompMat(WinPP, 0) Then
                SwapFlag = 0
            Else
                SwapFlag = 1
            End If
        Else 'NI recombinant
            NOPINI(0, SEventNumber) = 2: NOPINI(1, SEventNumber) = 1: NOPINI(2, SEventNumber) = 0
            DMiMa(INList(0)) = 2: DMiMa(INList(1)) = 1: DMiMa(INList(2)) = 0
            If INList(0) = CompMat(WinPP, 0) Then
                SwapFlag = 0
            Else
                SwapFlag = 1
            End If
        End If
    Else 'there has been no change in the relative tree positions of the three sequences
         NOPINI(0, SEventNumber) = WinPP
         DMiMa(WinPP) = 0
         If FAMat(WinPP, CompMat(WinPP, 0)) < FAMat(CompMat(WinPP, 0), CompMat(WinPP, 1)) Or FAMat(WinPP, CompMat(WinPP, 1)) < FAMat(CompMat(WinPP, 0), CompMat(WinPP, 1)) Then 'if winpp is inlyer
            If FAMat(WinPP, CompMat(WinPP, 0)) < FAMat(WinPP, CompMat(WinPP, 1)) Then
                'compmat 0 the other inlyer
                NOPINI(1, SEventNumber) = CompMat(WinPP, 1) 'minor parent the outlyer
                NOPINI(2, SEventNumber) = CompMat(WinPP, 0) 'major parent the other inlyer
                DMiMa(CompMat(WinPP, 1)) = 1
                DMiMa(CompMat(WinPP, 0)) = 2
                If SMat(WinPP, CompMat(WinPP, 0)) > 0 Or SMat(WinPP, CompMat(WinPP, 1)) > 0 Or SMat(CompMat(WinPP, 1), CompMat(WinPP, 0)) > 0 Then
                    If FMat(WinPP, CompMat(WinPP, 0)) > 0 Or FMat(WinPP, CompMat(WinPP, 1)) > 0 Or FMat(CompMat(WinPP, 1), CompMat(WinPP, 0)) > 0 Then

                        If FMat(WinPP, CompMat(WinPP, 0)) / (FMat(WinPP, CompMat(WinPP, 0)) + FMat(WinPP, CompMat(WinPP, 1)) + FMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) < SMat(WinPP, CompMat(WinPP, 0)) / (SMat(WinPP, CompMat(WinPP, 0)) + SMat(WinPP, CompMat(WinPP, 1)) + SMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) Then
                        'the sequences are further apart in the recombinant region
                            'minor parent actually unknown
                            MIF = 2
                        Else 'the sequences are closer together in the recombinant region
                            'major parent actually unknown
                            MIF = 1
                        End If
                    Else
                        MIF = 2
                    End If
                Else
                    MIF = 2
                End If
            Else
                'compmat 1 the other inlyer
                NOPINI(1, SEventNumber) = CompMat(WinPP, 0) 'minor parent the outlyer
                NOPINI(2, SEventNumber) = CompMat(WinPP, 1) 'major parent the other inlyer
                DMiMa(CompMat(WinPP, 1)) = 2
                DMiMa(CompMat(WinPP, 0)) = 1
                If SMat(WinPP, CompMat(WinPP, 0)) > 0 Or SMat(WinPP, CompMat(WinPP, 1)) > 0 Or SMat(CompMat(WinPP, 1), CompMat(WinPP, 0)) > 0 Then
                    If FMat(WinPP, CompMat(WinPP, 0)) > 0 Or FMat(WinPP, CompMat(WinPP, 1)) > 0 Or FMat(CompMat(WinPP, 1), CompMat(WinPP, 0)) > 0 Then
                        If FMat(WinPP, CompMat(WinPP, 1)) / (FMat(WinPP, CompMat(WinPP, 0)) + FMat(WinPP, CompMat(WinPP, 1)) + FMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) < SMat(WinPP, CompMat(WinPP, 1)) / (SMat(WinPP, CompMat(WinPP, 0)) + SMat(WinPP, CompMat(WinPP, 1)) + SMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) Then
                        'the sequences are further apart in the recombinant region
                            'minor parent actually unknown
                            MIF = 2
                        Else 'the sequences are closer together in the recombinant region
                            'major parent actually unknown
                            MIF = 1
                        End If
                    Else
                        MIF = 2
                    End If
                Else
                    MIF = 1
                End If
            End If
         Else 'if winpp is outlyer
            NOPINI(1, SEventNumber) = CompMat(WinPP, 0) 'minor parent one inlyer
            NOPINI(2, SEventNumber) = CompMat(WinPP, 1) 'major parent the other inlyer
            DMiMa(CompMat(WinPP, 1)) = 2
            DMiMa(CompMat(WinPP, 0)) = 1
            If (SMat(WinPP, CompMat(WinPP, 0)) + SMat(WinPP, CompMat(WinPP, 1)) + SMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) > 0 Then
            
                If ((FMat(WinPP, CompMat(WinPP, 0)) + FMat(WinPP, CompMat(WinPP, 1))) / 2) / (FMat(WinPP, CompMat(WinPP, 0)) + FMat(WinPP, CompMat(WinPP, 1)) + FMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) < ((SMat(WinPP, CompMat(WinPP, 0)) + SMat(WinPP, CompMat(WinPP, 1))) / 2) / (SMat(WinPP, CompMat(WinPP, 0)) + SMat(WinPP, CompMat(WinPP, 1)) + SMat(CompMat(WinPP, 1), CompMat(WinPP, 0))) Then
                    'outlyer further in the recombinant region than it is in bKg
                    'minor parent actually unknown
                    MIF = 2
                Else
                     'outlyer closer in the recombinant region than it is in bKg
                     'major parent actually unknown
                     MIF = 1
                End If
            Else 'all the sequences are identical in the recombinant region
                MIF = 1
               ' If FMat(WinPP, CompMat(WinPP, 0)) > FMat(WinPP, CompMat(WinPP, 1)) And FMat(WinPP, CompMat(WinPP, 0)) > FMat(CompMat(WinPP, 1), CompMat(WinPP, 0)) Then
               '     'winpp is outlyer
               '
               ' ElseIf FMat(CompMat(WinPP, 0), CompMat(WinPP, 1)) > FMat(WinPP, CompMat(WinPP, 1)) And FMat(CompMat(WinPP, 0), CompMat(WinPP, 1)) > FMat(WinPP, CompMat(WinPP, 0)) Then
               '     'CompMat(WinPP, 0) is outlyer
               '
               ' Else
               '     'CompMat(WinPP, 1) is outlyer
               '
               ' End If
            End If
         End If
         If NOPINI(1, SEventNumber) = CompMat(WinPP, 0) Then
            SwapFlag = 0
         Else
            SwapFlag = 1
         End If
        
    End If
    'Call MakeRange(0, WinPP, RNum(), RList, WinnerPos, StartPos, EndPos)
    'Call CBrother(NextNo, SeqNum(), StrainSeq(), ISeqs(), FMat(), FAMat(), SAMat(), NumSeqs, StartPos, EndPos)
    ReDim Breaks(1, RNum(WinPP))
    Call MakeBreaks(BPos, EPos, DMiMa(), MIF, SwapFlag, MissIDFlag, CompMat(), WinPP, ISeqs(), FMat(), SMat(), FAMat(), SAMat(), RecombNo(), MinPair(), Breaks(), BestXOList(), RNum(), TraceSub(), RList(), WinnerPos())
    Call MakeAlternative(BPos, EPos, DMiMa(), MIF, SwapFlag, MissIDFlag, CompMat(), CompMat(WinPP, 0), ISeqs(), FMat(), SMat(), FAMat(), SAMat(), MinPair(), tBXOListMa(), RNum(), TraceSub(), RList(), WinnerPosMa())
    Call MakeAlternative(BPos, EPos, DMiMa(), MIF, SwapFlag, MissIDFlag, CompMat(), CompMat(WinPP, 1), ISeqs(), FMat(), SMat(), FAMat(), SAMat(), MinPair(), tBXOListMi(), RNum(), TraceSub(), RList(), WinnerPosMi())
    X = X
    
    
    If SwapFlag = 0 Then
        Call CopyBXOList(TraceSub(), BestXOListMa(), BCurrentXoverMa(), tBXOListMa(), tBCurrentXoverMa())
        Call CopyBXOList(TraceSub(), BestXOListMi(), BCurrentXoverMi(), tBXOListMi(), tBCurrentXoverMi())
        Call CopyExtraHits(SEventNumber, ExtraHitsMi(), tEHitsMi())
        Call CopyExtraHits(SEventNumber, ExtraHitsMa(), tEHitsMa())
        
    Else
        Call CopyBXOList(TraceSub(), BestXOListMi(), BCurrentXoverMi(), tBXOListMa(), tBCurrentXoverMa())
        Call CopyBXOList(TraceSub(), BestXOListMa(), BCurrentXoverMa(), tBXOListMi(), tBCurrentXoverMi())
        Call CopyExtraHits(SEventNumber, ExtraHitsMi(), tEHitsMa())
        Call CopyExtraHits(SEventNumber, ExtraHitsMa(), tEHitsMi())
        
    End If
    
    ReDim Preserve DScores(20, 2, SEventNumber)
    For X = 0 To 20
            DScores(X, 0, SEventNumber) = TDScores(X, WinPP)
    Next X
    If SwapFlag = 0 Then
        For X = 0 To 20
            DScores(X, 2, SEventNumber) = TDScores(X, CompMat(WinPP, 1))
            DScores(X, 1, SEventNumber) = TDScores(X, CompMat(WinPP, 0))
        Next X
    Else
        For X = 0 To 20
            DScores(X, 1, SEventNumber) = TDScores(X, CompMat(WinPP, 1))
            DScores(X, 2, SEventNumber) = TDScores(X, CompMat(WinPP, 0))
        Next X
    End If
    '49.011 5k perms
    '0.250 5K perms - using makeeventseqs
    
    'find best event in rlist
    For X = 0 To RNum(WinPP)
        If RList(WinPP, X) = ISeqs(WinPP) Then Exit For
    Next X
    ActualE = X
    ReDim tDaught(0, Nextno)
    For X = 0 To RNum(WinPP)
        tDaught(0, RList(WinPP, X)) = 1
    Next X
    Call GetAge(SEventNumber, tDaught(), 0, ISeqs(WinPP), ISeqs(CompMat(WinPP, 1)), ISeqs(CompMat(WinPP, 0)), FMat(), SMat(), PermDiffs(), PermValid(), AgeEvent())
    Call UpdateAgeScore(TraceSub(ISeqs(WinPP)), AgeEvent(1, SEventNumber), SEventNumber, BPos, EPos, AgeScore(), EventScore())
    
    'XX = AgeScore(1000, ISeqs(WinPP))
    If AllowConflict = 0 And AbortFlag = 0 Then
        If SEventNumber >= UBound(RepeatCycles, 1) Then
            UB = UBound(RepeatCycles, 1)
            ReDim Preserve RepeatCycles(UB + 1000)
        End If
        RepeatCycles(SEventNumber) = RepeatCycles(SEventNumber) + 1
        
        
        Dim TestC() As Long
        Call TestConflict(TestC(), WinPP, AgeEvent(0, SEventNumber), Fail, BPos, EPos, TraceSub(), AgeEvent(), AgeScore(), EventScore(), ISeqs(), BestXOList(), BCurrentXOver())
        If Fail > 0 Then '2,1,3,2,1,3,2,4,3,3,3,3,3,3,3,3,3,3,3,3
            'oSEventNumber = SEventNumber
            If RepeatCycles(Fail) <= MaxRepeatCycles Then
                
                Call Rewind(TestC(), BCurrentXOver(), BestXOList())
                ReDim MissingData(Len(StrainSeq(0)), Nextno)
                oAbortFlag = AbortFlag
                
                Call BuildFirstXOList(1, SPX, AgeScore(), EventScore(), MinSeqSize, JumpFlag, MissingData(), TraceSub(), Nextno, StepNo, Steps(), ExtraHits(), ExtraHitsMa(), ExtraHitsMi(), NOPINI(), Eventnumber, SEventNumber, BestXOList(), BCurrentXOver(), XOverList(), CurrentXOver(), BestXOListMi(), BCurrentXoverMi(), BestXOListMa(), BCurrentXoverMa(), Daught(), MinorPar(), MajorPar())
                AbortFlag = oAbortFlag
                Call RemoveAccepts(BCurrentXOver(), BestXOList())
                'If AbortFlag = 0 Then
                    GoTo RestartX
                'End If
            Else
                X = X
                'ADD WARNING OF POSSIBLE CONFLICT IF YOURE JUST PUSHING ON
            '    'SEventNumber = oSEventNumber
            '    For X = SEventNumber + 1 To UBound(RepeatCycles, 1)
            '        RepeatCycles(X) = 0
            '    Next X
            End If
        Else
            X = X
        End If
        
            
        
    End If
    
    'begin section 10****************************************************************
    
    
      'erase bits from recombinants
     
   'For X = 0 To RNum(WinPP)
   '     Breaks(0, X) = BPos
   '     Breaks(1, X) = EPos
   '
   'Next X
    Dummy = ModSeqNumY(BPos, EPos, Len(StrainSeq(0)), WinPP, RNum(0), Breaks(0, 0), RList(0, 0), SeqNum(0, 0), tSeqNum(0, 0), MissingData(0, 0))
    
    'If MissingData(5450, 12) = 1 Then
    '    XX = ISeqs(0)
    '    XX = ISeqs(1)
    '    XX = ISeqs(2)
    '    XX = SEventNumber
    '    For X = 0 To RNum(1)
    '        XX = RList(1, X)
    '    Next X
    'End If
    
    ' 0.190
    
    ReDim DoPairs(Nextno, Nextno), TCurrentXOver(Nextno), TempDone(Nextno + RNum(WinPP) + 1, UBound(PXOList, 2))
   
    
    'place all evidence of non winners from pxolist in tempxolist
    oRecombNo(100) = 0
    '*****************************************************************************************
    'must set up code to redim preserve temxolist
    For X = 0 To Nextno
        For Y = 1 To PCurrentXOver(X)
            DA = PXOList(X, Y).Daughter
            Mi = PXOList(X, Y).MinorP
            Ma = PXOList(X, Y).MajorP
            
            
            WinPPY = MakePairs(Nextno, DA, Ma, Mi, WinPP, RNum(0), RList(0, 0), DoPairs(0, 0))
            If PXOList(X, Y).Probability <= LowestProb Then
                If WinPPY = RNum(WinPP) + 1 Then
                    If TCurrentXOver(TraceSub(DA)) <= TCurrentXOver(TraceSub(Mi)) And TCurrentXOver(TraceSub(DA)) <= TCurrentXOver(TraceSub(Ma)) Then
                        TCurrentXOver(X) = TCurrentXOver(X) + 1
                        If TCurrentXOver(X) <= XOSize Then
                            
                            TempXOList(X, TCurrentXOver(X)) = PXOList(X, Y)
                            If TCurrentXOver(X) <= UBound(TempDone, 2) And Y <= UBound(DoneSeq, 2) Then
                                TempDone(X, TCurrentXOver(X)) = DoneSeq(X, Y)
                            End If
                        Else
                            TCurrentXOver(X) = TCurrentXOver(X) - 1
                        End If
                    ElseIf TCurrentXOver(TraceSub(Mi)) <= TCurrentXOver(TraceSub(DA)) And TCurrentXOver(TraceSub(Mi)) <= TCurrentXOver(TraceSub(Ma)) Then
                        TCurrentXOver(TraceSub(Mi)) = TCurrentXOver(TraceSub(Mi)) + 1
                        If TCurrentXOver(TraceSub(Mi)) <= XOSize Then
                            TempXOList(Mi, TCurrentXOver(TraceSub(Mi))) = PXOList(X, Y)
                            TempXOList(Mi, TCurrentXOver(TraceSub(Mi))).Daughter = Mi
                            TempXOList(Mi, TCurrentXOver(TraceSub(Mi))).MinorP = DA
                            If TCurrentXOver(TraceSub(Mi)) <= UBound(TempDone, 2) And Y <= UBound(DoneSeq, 2) Then
                                TempDone(TraceSub(Mi), TCurrentXOver(TraceSub(Mi))) = DoneSeq(X, Y)
                            End If
                        Else
                            TCurrentXOver(TraceSub(Mi)) = TCurrentXOver(TraceSub(Mi)) - 1
                        End If
                        
                    Else
                        TCurrentXOver(TraceSub(Ma)) = TCurrentXOver(TraceSub(Ma)) + 1
                        If TCurrentXOver(TraceSub(Ma)) <= XOSize Then
                            TempXOList(Ma, TCurrentXOver(TraceSub(Ma))) = PXOList(X, Y)
                            TempXOList(Ma, TCurrentXOver(TraceSub(Ma))).Daughter = Ma
                            TempXOList(Ma, TCurrentXOver(TraceSub(Ma))).MajorP = DA
                            If TCurrentXOver(TraceSub(Ma)) <= UBound(TempDone, 2) And Y <= UBound(DoneSeq, 2) Then
                                TempDone(Ma, TCurrentXOver(TraceSub(Ma))) = DoneSeq(X, Y)
                            End If
                        Else
                            TCurrentXOver(TraceSub(Ma)) = TCurrentXOver(TraceSub(Ma)) - 1
                        End If
                    End If
                    oRecombNo(100) = oRecombNo(100) + 1
                    oRecombNo(PXOList(X, Y).ProgramFlag) = oRecombNo(PXOList(X, Y).ProgramFlag) + 1
                End If
            End If
        Next Y
    Next X
    
    For z = 0 To Nextno
        For X = 0 To RNum(WinPP)
            If DoPairs(RList(WinPP, X), z) = 1 Then
                For Y = 0 To RNum(WinPP)
                    DoPairs(RList(WinPP, Y), z) = 1
                    DoPairs(z, RList(WinPP, Y)) = 1
                Next Y
                Exit For
                'exit dot
            End If
        Next X
        
    Next z
    'If SEventNumber = 15 Then
    '    zz = 0
    '    For X = 0 To NextNo
    '        For Y = 0 To NextNo
    '            zz = zz + DoPairs(X, Y)
    '        Next Y
    '    Next X
    '    zz = zz
    '    For X = 0 To NextNo
    '        For Y = 1 To CurrentXover(X)
    '            If PXOList(X, Y).Daughter = X Then X = X
    '        Next Y
    '    Next X
    '    '158'162
    '    XX = ISeqs(WinPP)
    '    '25
    '
    'End If
    
    'erase all rescans from redolist
    If RedoListSize > 0 Then
        'mark removals
        'If X = X Then
            Dummy = MarkRemovals(Nextno, WinPP, RedoListSize, RedoList(0, 0), RNum(0), RList(0, 0), DoPairs(0, 0))
        'Else
        '    For X = 0 To Redolistsize
        '        For Y = 1 To 3
        '            For Z = 0 To RNum(WinPP)
        '                If Rlist(WinPP, Z) = RedoList(Y, X) Then
        '                    If Y = 1 Then
        '                        If DoPairs(RedoList(2, X), RedoList(3, X)) = 1 Then Exit For
        '                    ElseIf Y = 2 Then
        '                        If DoPairs(RedoList(1, X), RedoList(3, X)) = 1 Then Exit For
        '                    ElseIf Y = 3 Then
        '                         If DoPairs(RedoList(1, X), RedoList(2, X)) = 1 Then Exit For
        '                    End If
        '
        '                End If
        '            Next Z
        '            If Z <= RNum(WinPP) Then
        '                RedoList(0, X) = -1
        '                Exit For
        '            End If
        '        Next Y
        '    Next X
        'End If
        'clean up redolist
        Dummy = CleanRedoList(RedoListSize, RedoList(0, 0))
    End If
    
    
    'end section 10************************************************
    '9.844 5k perms
    '1.041
    '0.34 with makepairs
    'begin section 11************************************************
     
    '5.708 10 times more events
    '3.525 with makepairs
     
     'SSS = GetTickCount
    
    ' For i = 0 To 500
     'begin section 11.1************************************************
     
    
    ReDim CurrentXOver(Nextno), SLookUp(1, Nextno + 1), SLookUpNum(1)
    ReDim DonePVCO(AddNum - 1, Nextno), MaxXOP(AddNum - 1, Nextno)
    XX = oRecombNo(8)
    Call SignalCount(PXOList(), PCurrentXOver())
    Call UpdateRecNums(SEventNumber)
    'Redolistsize = 0
    For X = 0 To AddNum - 1
        For Y = 0 To Nextno
            DonePVCO(X, Y) = -1
        Next Y
    Next X
     
     If IndividualB > -1 Then
            
            Seq1 = ISeqs(WinPP)
            
            If MaskSeq(Seq1) = 0 Then
                SLookUpNum(1) = 1
                For X = 0 To Nextno
                    
                    If X <> Seq1 And MaskSeq(X) < 2 And ActualSeqSize(X) > MinSeqSize Then
                        If TraceSub(X) <> IndividualA And TraceSub(X) <> IndividualB Then
                            SLookUpNum(1) = SLookUpNum(1) + 1
                            SLookUp(1, SLookUpNum(1)) = X
                        Else
                            If X = IndividualA Or X = IndividualB Then
                                SLookUpNum(0) = 1
                                If TraceSub(Seq1) = IndividualA Then
                                    SLookUp(0, 1) = IndividualB
                                Else
                                    SLookUp(0, 1) = IndividualA
                                End If
                            End If
                        End If
                    End If
                Next X
            Else
                SLookUpNum(0) = 1
                SLookUp(0, 1) = IndividualA
                SLookUpNum(1) = 2
                SLookUp(1, 2) = IndividualB
            End If
        
     ElseIf IndividualA > -1 Then
        
            Seq1 = TraceSub(ISeqs(WinPP))
            If MaskSeq(Seq1) = 0 Then
                For X = 0 To Nextno
                    
                    If TraceSub(X) <> Seq1 And MaskSeq(X) < 2 And ActualSeqSize(X) > MinSeqSize Then
                        SLookUpNum(0) = SLookUpNum(0) + 1
                        SLookUp(0, SLookUpNum(0)) = X
                        SLookUpNum(1) = SLookUpNum(1) + 1
                        SLookUp(1, SLookUpNum(1)) = X
                    End If
                Next X
            Else
                For X = 0 To Nextno
                    If TraceSub(X) = IndividualA Then
                        SLookUpNum(0) = SLookUpNum(0) + 1
                        SLookUp(0, 1) = X
                    End If
                Next X
                
                SLookUpNum(1) = 1
                For X = 0 To Nextno
                    If X <> Seq1 And MaskSeq(X) < 2 And TraceSub(X) <> IndividualA And ActualSeqSize(X) > MinSeqSize Then
                        SLookUpNum(1) = SLookUpNum(1) + 1
                        SLookUp(1, SLookUpNum(1)) = X
                    End If
                Next X
            End If
     Else
        For X = 0 To Nextno
            If MaskSeq(X) = 0 Then
                SLookUpNum(0) = SLookUpNum(0) + 1
                SLookUp(0, SLookUpNum(0)) = X
                SLookUpNum(1) = SLookUpNum(1) + 1
                SLookUp(1, SLookUpNum(1)) = X
            End If
        Next X
     End If
     
     'end section 11.1 ************************************************
     '0.060 5K perms
     
     MCCorrectX = (RNum(WinPP) + 1) * (Nextno + 1) * (Nextno) / 2
     If SSOutlyerFlag = 2 Then
        Call GetOutie
        oSeq = Outie
     End If
     ReDim TSub(Nextno)
     For X = 0 To Nextno
        TSub(X) = X
     Next X
     B = 0
     'XX = StraiName(ISeqs(WinPP))
     Dim SAll As Long
     SAll = Abs(GetTickCount)
     For WinPPY = 0 To RNum(WinPP)
      'For Seq1 = 0 To Nextno
        'scan seqx against all the rest
        ' ie similar to individualA scan
        
        Seq1 = RList(WinPP, WinPPY)
        
        GoOn = 0
        If IndividualA = TraceSub(ISeqs(WinPP)) Or IndividualB = TraceSub(ISeqs(WinPP)) Then
            If TraceSub(Seq1) = IndividualA Or TraceSub(Seq1) = IndividualB Then
                GoOn = 1
            End If
        Else
            GoOn = 1
        End If
        
        GoOn = 1
        
        If GoOn = 1 And ActualSeqSize(Seq1) > MinSeqSize And MaskSeq(Seq1) = 0 Then
            'For Seq2 = Seq1 + 1 To Nextno
            For G = 1 To SLookUpNum(0)
                Seq2 = SLookUp(0, G)
                If ActualSeqSize(Seq2) > MinSeqSize Then
                    For X = 0 To WinPPY
                        If Seq2 = RList(WinPP, X) Then Exit For
                    Next X
                
                    'X = WinPPY + 1
                    If X = WinPPY + 1 Then
                        For H = G + 1 To SLookUpNum(1)
                            Seq3 = SLookUp(1, H)
                            If ActualSeqSize(Seq3) > MinSeqSize Then
                                For X = 0 To WinPPY
                                    If Seq3 = RList(WinPP, X) Then Exit For
                                Next X
                                If X = WinPPY + 1 Then
                                    If DoPairs(Seq2, Seq3) = 1 Then '
                                        If PermValid(Seq1, Seq2) > MinSeqSize And PermValid(Seq1, Seq3) > MinSeqSize And PermValid(Seq2, Seq3) > MinSeqSize Then
                                            'If X = X Or SubValid(Seq1, Seq2) > 20 And SubValid(Seq1, Seq3) > 20 And SubValid(Seq2, Seq3) > 20 Then
    
                                                If TraceSub(Seq1) <> TraceSub(Seq2) And TraceSub(Seq1) <> TraceSub(Seq3) And TraceSub(Seq2) <> TraceSub(Seq3) Then
                                                    
                                                    If DoScans(0, 0) = 1 Then Call XOver(SeqNum(), Seq1, Seq2, Seq3, 0)
                                                        
                                                    If DoScans(0, 1) = 1 Then Call GCXoverD(0)
                                                    If DoScans(0, 2) = 1 Then
                                                            'BSStepsize = BSStepsize
                                                        Call BSXoverS(BackUpNextno, Seq1, Seq2, Seq3, TraceSub(), MissingData())
                                                    End If
                                                     
                                                    If DoScans(0, 3) = 1 Then Call MCXoverF(0, 0, 0)
                                                    X = X
                                                    If DoScans(0, 4) = 1 Then
                                                        tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                                                                
                                                        Call CXoverA(0, 0, 0)
                                                                
                                                        Seq1 = tseq2: Seq2 = tseq3: Seq3 = tseq1
                                                                
                                                        Call CXoverA(0, 0, 0)
                                                                
                                                        Seq1 = tseq3: Seq2 = tseq1: Seq3 = tseq2
                                                                
                                                        Call CXoverA(0, 0, 0)
                                                                
                                                        Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                                                    End If
                                                    If DoScans(0, 8) = 1 Then
                                                        tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                                                                
                                                        Call TSXOver(0)
                                                                
                                                        Seq1 = tseq2: Seq2 = tseq3: Seq3 = tseq1
                                                                
                                                        Call TSXOver(0)
                                                                
                                                        Seq1 = tseq3: Seq2 = tseq1: Seq3 = tseq2
                                                                
                                                        Call TSXOver(0)
                                                                
                                                        Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                                                    End If
                                                    If DoScans(0, 5) = 1 Then
                                                        oSeq1 = Seq1: oSeq2 = Seq2: oSeq3 = Seq3
                                                        Call SSXoverC(0, WinNum, SeqMap(), ZPScoreHolder(), ZSScoreHolder(), CorrectP, oSeq, PermSScores(), PermPScores(), SScoreHolder(), PScoreHolder(), TraceSub(), SeqScore3(), MeanPScore(), SDPScore(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP())
                                                        Seq1 = oSeq1: Seq2 = oSeq2: Seq3 = oSeq3
                                                    End If
                                                    B = B + 1
                                                    
                                                End If
                                            'End If
                                        End If
                                        
                                    'Else
                                    '    X = X
                                    End If
                                End If
                            End If
                            ET = Abs(GetTickCount)
                            If Abs(ET - LT) > 500 Then
                                LT = ET
                                Form1.SSPanel1.Caption = Trim(Str(B)) & " of " & Trim(Str(MCCorrectX)) & " triplets reexamined"
                                If TT - elt > 2000 Then
                                            elt = ET
                                            If oTotRecs > 0 Then
                                                pbv = (1 - (oRecombNo(100) ^ 0.4 / oTotRecs ^ 0.4)) * 100
                                                If pbv > Form1.ProgressBar1 Then
                                                    Form1.ProgressBar1 = pbv
                                                End If
                                            End If
                                            
                                End If
                                
                                
                                
                                
                                If AbortFlag = 1 Then
                                    WinPPY = Nextno
                                    G = Nextno
                                    H = Nextno
                                End If
                                UpdateRecNums (SEventNumber)
                                Form1.Label57(0).Caption = DoTimeII(Abs(ET - STime))
                                Form1.Label50(0).Caption = DoTimeII(Abs(ET - SAll) * Abs(TimeFract(0)) + MethodTime(0))
                                Form1.Label51(0).Caption = DoTimeII(Abs(ET - SAll) * Abs(TimeFract(1)) + MethodTime(1))
                                Form1.Label66(0).Caption = DoTimeII(Abs(ET - SAll) * Abs(TimeFract(2)) + MethodTime(2))
                                Form1.Label67(0).Caption = DoTimeII(Abs(ET - SAll) * Abs(TimeFract(3)) + MethodTime(3))
                                Form1.Label68(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(4) + MethodTime(4))
                                Form1.Label69(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(5) + MethodTime(5))
                                Form1.Label5(1).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(8) + MethodTime(8))
                            End If
                        Next H
                    End If
                End If
            Next G
        End If
    Next WinPPY
   
    ET = Abs(GetTickCount)
    For X = 0 To AddNum - 1
        MethodTime(X) = Abs(MethodTime(X))
    Next X
    Form1.Label50(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(0) + MethodTime(0))
    Form1.Label51(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(1) + MethodTime(1))
    Form1.Label66(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(2) + MethodTime(2))
    Form1.Label67(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(3) + MethodTime(3))
    Form1.Label68(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(4) + MethodTime(4))
    Form1.Label69(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(5) + MethodTime(5))
    Form1.Label5(1).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(8) + MethodTime(8))
    For X = 0 To AddNum - 1
        MethodTime(X) = MethodTime(X) + Abs(ET - SAll) * TimeFract(X)
    Next X
    'If doscans(0,5) = 1 Then
    '    Call SSXoverD(0, WinPP, MinSeqSize, ActualSeqSize(), RList(), RNum(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP(), TraceSub(), SimSeqNum())
    ''    SSOutlyerFlag = pSSOutlyerFlag
    'End If
    'Form1.SSPanel1.Caption = Trim(Str(mccorrectx)) & " of " & Trim(Str(mccorrectx)) & " triplets reexamined"
    UpdateRecNums (SEventNumber)
    
   'Next i
 'ee = GetTickCount
  '  tt = ee - SSS
  '  X = X
    
    'maxchi
    '14.46 500 perms
    '7.860 (dsktop)
    '5.828 with growmchiwin
    '5.172 - not bothering with twins
    '5.878 (laptop) - rearranging xolist less
    '2.874 (laptop) - stopping when mpv doesn't make the cut - there may be a screwup
    '6.269 - more events
    '5.799 - destroy peaks
    '5.168 - findside
    '5.138 - != findsubseq
    '5.017 improvements in findsubseqc
    '4.837 improvements in winscorecalc
    '4.637 - taking out a criticaldiff test in winscorecalc.
    '4.376 - using criticalval in chipvals
    '10.212 with looking for overlapping missing data
    '9.835 only scrubbing negative peaks (ie not the whole region)
    '4.126 using banwins
    '3.906 using destroypeak
    '3.805 more efficient calculation chivals in growmchiwin
    '3.710
    '3.495 with hiseqs
    
    'X = X
  '  '24.986 500 perms
    '12.378 500 perms (using highenough)
    '11.647
    '11.186
    '14.801 with splits
    '11.707 with makesubprob
    'end section 11********************************************
    'GENECONV results
    '41.375 5K perms ignor indels (desktop)
    '55.109 (laptop)
    '53.397 (laptop) better array dimentioning
    '250.801 -use indels (laptop)
    
    
    '8.250 1k perms
    
    
    '326.439
    '115.005 screening only "relevant" pairs
    '86.704 - Improvement to findsubseq and rearrangement of ifs etc in the main loop
    '82.672 - inclusion of corrcetprob in probcalc
    '62.219 - improvements in xover
    '61.248 - taking out pointless cycles
    '60.938 - not using pointers in findnext
    
    '19.907 1K perms (0.282 without callinq xover)
    '19.094 1K perms - slight improvement in findsubseq
    '17.438 1K perms - rearrangement of ifs etc in main loop.
    '16.406 1K perms - inclusion of corrcetprob in probcalc
    '14.594 1K perms - only accessing xoverlist at the end after p-val is confirmed
    '12.859 1K perms - xoverwin and lenstrainseq calculated outside loop
    
    
    'Begin section 12*****************************************
    
    'Replace the missing bits
    
    oRnum(WinPP) = RNum(WinPP)
    ReDim oRlist(RNum(WinPP)), oBreaks(1, RNum(WinPP))
    For X = 0 To RNum(WinPP)
        oRlist(X) = RList(WinPP, X)
        oBreaks(0, X) = Breaks(0, X)
        oBreaks(1, X) = Breaks(1, X)
    Next X
    Dummy = RebuildSeqNum(Len(StrainSeq(0)), WinPP, RNum(0), Breaks(0, 0), RList(0, 0), SeqNum(0, 0), tSeqNum(0, 0))
    
    
    'End section 12*********************************************
    '9.163 5K perms
    '0.010
    
    'Befor removing any events from RList it is very important to
    'first take note of which distances do and do not need to be recalculated.
    'distances between any sequence and all sequences in rlist at this point must
    'all be recalculated.
    
    'begin section 13*********************************************
    
    Dummy = MakeUninvolved(WinPP, Nextno, UnInvolved(0), RNum(0), RList(0, 0))
    
    'remove events that were not actually detected
    Dummy = StripUnfound2(WinPP, AddNum, RNum(0), WinnerPos(0, 0), RList(0, 0), Breaks(0, 0))
    
     
    'end section 13*********************************************
    '0.110
    '0.007
     
    'begin section 14*********************************************
    
    'ReDim TCurrentXOver(Nextno)
    'place evidence in xoverlist into tempxolist
    For X = 0 To Nextno
           For Y = 1 To CurrentXOver(X)
                TCurrentXOver(X) = TCurrentXOver(X) + 1
                If TCurrentXOver(X) > UBound(TempXOList, 2) Then
                    XOSize = TCurrentXOver(X) + 10
                    On Error Resume Next
                    TXOS = UBound(TempXOList, 1)
                    On Error GoTo 0
                    ReDim Preserve TempXOList(TXOS, XOSize)
                    UB = UBound(TempDone, 1)
                    ReDim Preserve TempDone(UB, XOSize)
                End If
                TempXOList(X, TCurrentXOver(X)) = XOverList(X, Y)
                 
                                        
           Next Y
    Next X
    
    
    
    
    MXOSize = 0
    For X = 0 To Nextno
        If TCurrentXOver(X) > MXOSize Then
        MXOSize = TCurrentXOver(X)
        End If
    Next X
    XOSize = MXOSize
    'end section 14*********************************************
    '0.140
    
    '0.721 10 times more events
    'begin section 15*********************************************
    
    'now redim permxover and copy tempxover to permxover
     
    ReDim PXOList(Nextno + RNum(WinPP) + 1, XOSize), PCurrentXOver(Nextno + RNum(WinPP) + 1), NumRecsI(Nextno)
    
    For X = 0 To Nextno
        PCurrentXOver(X) = TCurrentXOver(X)
        For Y = 1 To TCurrentXOver(X)
            
            PXOList(X, Y) = TempXOList(X, Y)
           
            DA = PXOList(X, Y).Daughter
            Ma = PXOList(X, Y).MajorP
            Mi = PXOList(X, Y).MinorP
            NumRecsI(TraceSub(DA)) = NumRecsI(TraceSub(DA)) + 1
            NumRecsI(TraceSub(Ma)) = NumRecsI(TraceSub(Ma)) + 1
            NumRecsI(TraceSub(Mi)) = NumRecsI(TraceSub(Mi)) + 1
        Next Y
    Next X
  
    'end section 15*********************************************
    '0.531 5K perms
    '0.906
    
    '3.258 10X more events
    'remove recombinant region from winner and make extra sequences
    
    'Find out whether it is necessary to redo former events in light of the fact that
    'the current event was previously used as a parent in the recombinant region.
    
    
    
    ReDim Preserve MissingData(Len(StrainSeq(0)), Nextno)
    xNextno = Nextno
    
    Nextno = Nextno + RNum(WinPP) + 1
    
    
    
    'begin section 16*********************************************
    
    ReDim Preserve GrpMaskSeq(Nextno), TraceSub(Nextno), NumRecsI(Nextno), SeqNum(Len(StrainSeq(0)), Nextno), UsedPar(Len(StrainSeq(0)), Nextno)
    ReDim Preserve SimSeqNum(Len(StrainSeq(0)), Nextno), SubMaskSeq(Nextno), ActualSeqSize(Nextno), MissingData(Len(StrainSeq(0)), Nextno)
    
    If Nextno > UBound(Relevant2, 2) Then
        ReDim Preserve Relevant2(2, Nextno)
    End If
    
    
    For X = 0 To RNum(WinPP)
        TraceSub(Nextno - RNum(WinPP) + X) = TraceSub(RList(WinPP, X))
    Next X
    
    'Make extra sequences and delete bits from recombinants
    Dummy = ModSN(Nextno, Len(StrainSeq(0)), BPos, EPos, WinPP, RNum(0), RList(0, 0), Breaks(0, 0), SeqNum(0, 0), MissingData(0, 0))
    Dummy = ModSeqNumZ(Nextno, Len(StrainSeq(0)), BPos, EPos, WinPP, oRnum(0), oRlist(0), oBreaks(0, 0), SeqNum(0, 0), MissingData(0, 0))
    
   
    'add to steps
    For X = 0 To RNum(WinPP)
        GoOn = 0
        For A = 0 To AddNum - 1
            If WinnerPos(X, A) > 0 Then
                'If StepNo = 995 Then
                '    X = X
                'End If
                Steps(0, StepNo) = 1 'ie create a new sequence ....
                Steps(1, StepNo) = RList(WinPP, X) 'using this seqence.....
                Steps(2, StepNo) = Breaks(0, X) 'from this position....
                Steps(3, StepNo) = Breaks(1, X) 'to this position....
                Steps(4, StepNo) = SEventNumber + 1
                StepNo = StepNo + 1
                
                UB = UBound(Steps, 2)
                If StepNo > UB Then
                    ReDim Preserve Steps(4, UB + 100)
                End If
                GoOn = 1
                Exit For
            End If
        Next A
        If GoOn = 0 Or X = 12345 Then
            Steps(0, StepNo) = 1 'ie create a new sequence ....
            Steps(1, StepNo) = RList(WinPP, X) 'using this seqence.....
            Steps(2, StepNo) = Breaks(0, X) 'from this position....
            Steps(3, StepNo) = Breaks(1, X) 'to this position....
            Steps(4, StepNo) = SEventNumber + 1
            StepNo = StepNo + 1
            
            UB = UBound(Steps, 2)
            If StepNo > UB Then
                ReDim Preserve Steps(4, UB + 100)
            End If
        End If
    Next X
    For X = 0 To oRnum(WinPP)
        'If StepNo = 995 Then
        '            X = X
        '        End If
        'XX = StraiName(ISeqs(WinPP))
        Steps(0, StepNo) = 2 'ie delete a bit of sequence ....
        Steps(1, StepNo) = oRlist(X)  'from this seqence.....
        Steps(2, StepNo) = oBreaks(0, X) 'from this position....
        Steps(3, StepNo) = oBreaks(1, X) 'to this position....
        Steps(4, StepNo) = SEventNumber + 1
        StepNo = StepNo + 1
        
        UB = UBound(Steps, 2)
        If StepNo > UB Then
            ReDim Preserve Steps(4, UB + 100)
        End If
    Next X
    
   
   
    'end section 16*********************************************
    '8.953 5K perms
    '0.090
    
    'Update Actualseqsize
    'begin section 17*********************************************
    
    Dummy = MakeActualSeqSize(Len(StrainSeq(0)), Nextno, WinPP, RNum(0), RList(0, 0), ActualSeqSize(0), SeqNum(0, 0))
    'If X = X And SEventNumber = 70 Then
    '    Open "seqsize.csv" For Output As #1
    '    For X = 0 To PermNextNo
    '    'For X = PermNextNo To 0 Step -1
    '        Print #1, ActualSeqSize(X)
    '    Next X
    '    Close #1
    '    X = X
    'End If
    'end section 17*********************************************
    '25.486 5K perms
    '12.228 better array use
    '0.190 - using make actualseqsize
    
    'Make sure that
    '(1) ISeqs(rlist(winpp,winppy)) and nextno are not below the minimum acceptable
    'size setting - if it does don't add it
    '(2) nextno does not exceed the maxseqnum setting
    '    If it does boot out the smallest current seq
    'GoOn = 1
    
    'begin section 18*********************************************
    
    If BackUpNextno < Nextno Then
        ReDim Preserve SeqNum(Len(StrainSeq(0)), Nextno), UsedPar(Len(StrainSeq(0)), Nextno), SimSeqNum(Len(StrainSeq(0)), Nextno), TraceSub(Nextno), StraiName(Nextno), MaskSeq(Nextno)
        'Update the treedistance used for RDP and SiScan
        ReDim TreeDistance(Nextno, Nextno)
        ReDim Distance(Nextno, Nextno)
    End If
   
   
    '0.040
    '1.828 500K reps
    
    'Update distances
    'we have values already calculated for:
    '(1) a set fo sequences that are uninvolved in recombination - ie sequences outside rlist
    'and sequences =< xnextno
    
    ReDim TDiffs(Nextno, Nextno), TValid(Nextno, Nextno)
     
    Dummy = MakeTDiffs(Nextno, xNextno, TDiffs(0, 0), TValid(0, 0), PermDiffs(0, 0), PermValid(0, 0), UnInvolved(0), ActualSeqSize(0))
    'XX = UnInvolved(9)
    'XX = TValid(9, 106)
    'Why does his take so long? - it could be quicker (but more complicated) to
    'use sebdiffs and subvalids for sequences > last nexto with
    'sequences < lastnextno
    'I could also speed things up by keeping trackof te first and
    'last characters in the sequence strings and only looking between these
    
    ReDim Preserve UnInvolved(Nextno)
    ReDim PermDiffs(Nextno, Nextno), PermValid(Nextno, Nextno)
    
    
    For X = xNextno + 1 To Nextno
        XX = RNum(WinPP)
        UnInvolved(X) = 0
    Next X
    
    'ss = GetTickCount
    'For X = 0 To 5000
    AvDst = 0
    UDst = DistanceCalcW(Nextno, Len(StrainSeq(0)) + 1, TDiffs(0, 0), TValid(0, 0), SeqNum(0, 0), Distance(0, 0), AvDst, UnInvolved(0))
    DistanceFlag = 1
    
    
    'Next X
    'ee = GetTickCount
    'tt = ee - ss
    'X = X
    '6.549
    
    Dummy = MakePermDiffs(Nextno, MinSeqSize, TDiffs(0, 0), TValid(0, 0), PermDiffs(0, 0), PermValid(0, 0), Distance(0, 0))
    
    
    'Call CheckDists(SeqNum(), PermValid())
    'If NextNo > 105 Then
    '    If PermValid(9, 106) = 0 Then
    '    XX = NextNo
    '    End If
    'End If
   ' XX = 0
   ' For Z = 1 To Len(StrainSeq(0))
   '     If SeqNum(Z, 9) <> 46 Then
   '         If SeqNum(Z, 106) <> 46 Then
   '
   '             XX = XX + 1
   '         End If
   '     End If
   '
   ' Next Z
   ' X = X
    'end section 18*********************************************
    '13.189
    '7.010 using distancecacw
    '6.870
    '7.291
    '7.011 - using maketdist
    '6.589 - uing makepermvalid
    '5.703 -dsktop
    
   
    'begin section 19*********************************************
  '  ss = GetTickCount
  '  For X = 0 To 5000
    TreeDistFlag = 0
    Call UPGMA(0, 1)
    
  ' Next X
  ' ee = GetTickCount
  ' tt = ee - ss
  
    'end section 19*********************************************
    '2.824 5k perms
    '2.714 5 K perms
    '1.422 5K perms with addseqstoUPGMA
    '0.781 5k perms with treedist2
    'Begin section 20*********************************************
    
    'if the recombinant region fragment is added then scan it (ie seq1)
    'against the other sequences.
    ReDim DonePVCO(AddNum - 1, Nextno), MaxXOP(AddNum - 1, Nextno)
    ReDim CurrentXOver(Nextno)
    Call SignalCount(PXOList(), PCurrentXOver())
    Call UpdateRecNums(SEventNumber)
    For X = 0 To AddNum - 1
        For Y = 0 To Nextno
            DonePVCO(X, Y) = -1
        Next Y
    Next X
    
    'change this - xosize can get VERY VERY BIG- must work on ways to keep xosize as small
    'as possible but must also reset the xoverlist size to something smaller
    oxosize = XOSize
    XOSize = XOverListSize
    If Nextno > SNextno Then ReDim XOverList(Nextno, XOSize)
    XOSize = oxosize
    For X = Nextno - RNum(WinPP) To Nextno
        If TraceSub(X) <= Nextno Then
            MaskSeq(X) = MaskSeq(TraceSub(X)): GrpMaskSeq(X) = GrpMaskSeq(TraceSub(X))
        Else
            TraceSub(X) = X
        End If
    Next X
     
    For X = Nextno - RNum(WinPP) To Nextno
        If ActualSeqSize(X) > MinSeqSize Then Exit For
    Next X
        
        If X <= Nextno Then
            ReDim SLookUp(1, Nextno + 1), SLookUpNum(1)
            
            
            If IndividualB > -1 Then
                'Seq1 = Nextno - RNum(WinPP) + WinPPY
                For X = Nextno - RNum(WinPP) To Nextno
                    If TraceSub(X) = IndividualA Or TraceSub(X) = IndividualB Then Exit For
                Next X
                If X <= Nextno Then
                    Seq1 = TraceSub(X)
                    If Seq1 = IndividualA Then
                        SLookUp(0, 1) = IndividualB
                    Else
                        SLookUp(0, 1) = IndividualA
                    End If
                    SLookUpNum(1) = 1
                    For X = 0 To Nextno - RNum(WinPP) - 1
                        If TraceSub(X) <> Seq1 And MaskSeq(X) < 2 And ActualSeqSize(X) > MinSeqSize Then
                            If TraceSub(X) <> IndividualA And TraceSub(X) <> IndividualB Then
                                SLookUpNum(1) = SLookUpNum(1) + 1
                                SLookUp(1, SLookUpNum(1)) = X
                            End If
                        End If
                    Next X
                Else
                    SLookUpNum(0) = 1
                    SLookUp(0, 1) = IndividualA
                    SLookUpNum(1) = 2
                    SLookUp(1, 2) = IndividualB
                End If
               
            ElseIf IndividualA > -1 Then
                For X = Nextno - RNum(WinPP) To Nextno
                    If TraceSub(X) = IndividualA Then Exit For
                Next X
                If X <= Nextno Then
                    Seq1 = TraceSub(X)
                    For X = 0 To Nextno - RNum(WinPP) - 1
                        If TraceSub(X) <> Seq1 And MaskSeq(X) < 2 And ActualSeqSize(X) > MinSeqSize Then
                            SLookUpNum(0) = SLookUpNum(0) + 1
                            SLookUp(0, SLookUpNum(0)) = X
                            SLookUpNum(1) = SLookUpNum(1) + 1
                            SLookUp(1, SLookUpNum(1)) = X
                        End If
                    Next X
                Else
                    SLookUpNum(1) = 1
                    Seq1 = TraceSub(ISeqs(WinPP))
                    For X = 0 To Nextno - RNum(WinPP) - 1
                        If TraceSub(X) = IndividualA Then
                            SLookUpNum(0) = SLookUpNum(0) + 1
                            SLookUp(0, 1) = X
                        ElseIf TraceSub(X) <> Seq1 And MaskSeq(X) < 2 And ActualSeqSize(X) > MinSeqSize Then
                            SLookUpNum(1) = SLookUpNum(1) + 1
                            SLookUp(1, SLookUpNum(1)) = X
                        End If
                        
                    Next X
                       
                End If
                   
                       
                   
            Else
               'For WinPPY = 0 To RNum(WinPP)
                   For X = 0 To Nextno
                       If MaskSeq(X) = 0 Then
                           SLookUpNum(0) = SLookUpNum(0) + 1
                           SLookUp(0, SLookUpNum(0)) = X
                           SLookUpNum(1) = SLookUpNum(1) + 1
                           SLookUp(1, SLookUpNum(1)) = X
                       End If
                   Next X
               'Next WinPPY
            End If
        
            If SSOutlyerFlag = 2 Then
               Call GetOutie
               oSeq = Outie
            End If
            ReDim TSub(Nextno)
            For X = 0 To Nextno
               TSub(X) = X
            Next X
            
            B = 0
            MCCorrectX = (RNum(WinPP) + 1) * (Nextno + 1) * (Nextno) / 2
            SAll = GetTickCount
            For WinPPY = 0 To RNum(WinPP)
                Seq1 = Nextno - RNum(WinPP) + WinPPY
                GoOn = 0
                If IndividualA = TraceSub(ISeqs(WinPP)) Or IndividualB = TraceSub(ISeqs(WinPP)) Then
                    If TraceSub(Seq1) = IndividualA Or TraceSub(Seq1) = IndividualB Then
                        GoOn = 1
                    End If
                Else
                    GoOn = 1
                End If
                If GoOn = 1 And ActualSeqSize(Seq1) > MinSeqSize Then
                    For G = 1 To SLookUpNum(0)
                        Seq2 = SLookUp(0, G)
                        If ActualSeqSize(Seq2) > MinSeqSize Then
                            If Seq2 > SNextno Then
                                tseq2 = TraceSub(Seq2)
                            Else
                                tseq2 = Seq2
                            End If
                            
                            If SubValid(TraceSub(Seq1), tseq2) > 20 Then
                                For X = 0 To WinPPY
                                    If Seq2 = Nextno - RNum(WinPP) + X Then Exit For
                                Next X
                                If X = WinPPY + 1 Then
                                    
                                    For H = G + 1 To SLookUpNum(1)
                                        Seq3 = SLookUp(1, H)
                                        If ActualSeqSize(Seq3) > MinSeqSize Then
                                            For X = 0 To WinPPY
                                                If Seq3 = Nextno - RNum(WinPP) + X Then Exit For
                                            Next X
                                            
                                            
                                            
                                            If X = WinPPY + 1 Then
                                                If Seq3 > SNextno Then
                                                    tseq3 = TraceSub(Seq3)
                                                Else
                                                    tseq3 = Seq3
                                                End If
                                                If Seq2 > SNextno Then
                                                    tseq2 = TraceSub(Seq2)
                                                Else
                                                    tseq2 = Seq2
                                                End If
                                                If SubValid(tseq2, tseq3) > 20 And SubValid(TraceSub(Seq1), tseq3) > 20 Then
                                                    If TraceSub(Seq1) <> TraceSub(Seq2) And TraceSub(Seq1) <> TraceSub(Seq3) And TraceSub(Seq2) <> TraceSub(Seq3) Then
                                                        'XX = ActualSeqSize(Seq1)
                                                        'XX = ActualSeqSize(Seq2)
                                                        'XX = ActualSeqSize(Seq3)
                                                        'If ActualSeqSize(Seq1) < MinSeqSize Or ActualSeqSize(Seq2) < MinSeqSize Or ActualSeqSize(Seq1) < MinSeqSize Then
                                                        '    X = X
                                                        'End If
                                                        
                                                        If DoScans(0, 0) = 1 Then Call XOver(SeqNum(), Seq1, Seq2, Seq3, 0)
                                                        If DoScans(0, 1) = 1 Then Call GCXoverD(0)
                                                        
                                                        If DoScans(0, 2) = 1 Then
                                                            Call BSXoverS(BackUpNextno, Seq1, Seq2, Seq3, TraceSub(), MissingData())
                                                        End If
                                                        
                                                        If DoScans(0, 3) = 1 Then Call MCXoverF(0, 0, 0)
                                                       X = X
                                                        If DoScans(0, 4) = 1 Then
                                                            tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                                                                        
                                                            Call CXoverA(0, 0, 0)
                                                                        
                                                            Seq1 = tseq2: Seq2 = tseq3: Seq3 = tseq1
                                                                        
                                                            Call CXoverA(0, 0, 0)
                                                                        
                                                            Seq1 = tseq3: Seq2 = tseq1: Seq3 = tseq2
                                                                        
                                                            Call CXoverA(0, 0, 0)
                                                                        
                                                            Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                                                        End If
                                                        If DoScans(0, 8) = 1 Then
                                                            tseq1 = Seq1: tseq2 = Seq2: tseq3 = Seq3
                                                                    
                                                            Call TSXOver(0)
                                                                    
                                                            Seq1 = tseq2: Seq2 = tseq3: Seq3 = tseq1
                                                                    
                                                            Call TSXOver(0)
                                                                    
                                                            Seq1 = tseq3: Seq2 = tseq1: Seq3 = tseq2
                                                                    
                                                            Call TSXOver(0)
                                                                    
                                                            Seq1 = tseq1: Seq2 = tseq2: Seq3 = tseq3
                                                        End If
                                                        If DoScans(0, 5) = 1 Then
                                                            oSeq1 = Seq1: oSeq2 = Seq2: oSeq3 = Seq3
                                                            Call SSXoverC(0, WinNum, SeqMap(), ZPScoreHolder(), ZSScoreHolder(), CorrectP, oSeq, PermSScores(), PermPScores(), SScoreHolder(), PScoreHolder(), TraceSub(), SeqScore3(), MeanPScore(), SDPScore(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP())
                                                            Seq1 = oSeq1: Seq2 = oSeq2: Seq3 = oSeq3
                                                        End If
                                                    End If
                                                
                                                End If
                                                B = B + 1
                                            End If
                                                
                                            ET = Abs(GetTickCount)
                                            If Abs(ET - LT) > 500 Then
                                                
                                                LT = ET
                                                
                                                Form1.SSPanel1.Caption = Trim(Str(B)) & " of " & Trim(Str(MCCorrectX)) & " triplets reexamined"
                                                If Abs(ET - elt) > 2000 Then
                                                    elt = ET
                                                    If oTotRecs > 0 Then
                                                        pbv = (1 - (oRecombNo(100) ^ 0.4 / oTotRecs ^ 0.4)) * 100
                                                        If pbv > Form1.ProgressBar1 Then
                                                            Form1.ProgressBar1 = pbv
                                                        End If
                                                    End If
                                                    
                                                End If
                                                SAll = Abs(SAll)
                                                UpdateRecNums (SEventNumber)
                                                Form1.Label57(0).Caption = DoTimeII(Abs(ET - STime))
                                                Form1.Label50(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(0) + MethodTime(0))
                                                Form1.Label51(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(1) + MethodTime(1))
                                                Form1.Label66(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(2) + MethodTime(2))
                                                Form1.Label67(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(3) + MethodTime(3))
                                                Form1.Label68(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(4) + MethodTime(4))
                                                Form1.Label69(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(5) + MethodTime(5))
                                                Form1.Label5(1).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(8) + MethodTime(8))
                                                If AbortFlag = 1 Then
                                                    WinPPY = Nextno
                                                    G = Nextno
                                                    H = Nextno
                                                End If
                                            End If
                                        End If
                                    Next H
                                
                                End If
                            
                            End If
                        End If
                    Next G
                End If
            Next WinPPY
            SAll = Abs(SAll)
              
            
            ET = Abs(GetTickCount)
            Form1.Label50(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(0) + MethodTime(0))
            
            Form1.Label51(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(1) + MethodTime(1))
            Form1.Label66(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(2) + MethodTime(2))
            Form1.Label67(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(3) + MethodTime(3))
            Form1.Label68(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(4) + MethodTime(4))
            Form1.Label69(0).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(5) + MethodTime(5))
            Form1.Label5(1).Caption = DoTimeII(Abs(ET - SAll) * TimeFract(8) + MethodTime(8))
            For X = 0 To AddNum - 1
                MethodTime(X) = MethodTime(X) + Abs(ET - SAll) * TimeFract(X)
            Next X
            'If doscans(0,5) = 1 Then
            '    Call SSXoverD(1, WinPP, MinSeqSize, ActualSeqSize(), RList(), RNum(), Seq34Conv(), VRandConv(), VRandTemplate(), HRandTemplate(), TakenPos(), DG1(), DG2(), DoGroupS(), DoGroupP(), TraceSub(), SimSeqNum())
            '    SSOutlyerFlag = pSSOutlyerFlag
            'End If
            
            '194.199 - without any actual analysis
            '691.208 - rdp only
            '605.250 - maxchi only
            '1455.182 - geneconv alone
                '715.979 - findsubseqgca
                '1224.401 - getfrags
                '1384.200 - getmaxfrags
                '1455.052 - gccalcpval
                
            '2452.977 - with all
            
            'eee = GetTickCount
            'ttt = eee - SSS
            Form1.SSPanel1.Caption = Trim(Str(MCCorrectX)) & " of " & Trim(Str(MCCorrectX)) & " triplets reexamined"
            
            'Add events from xoverlist to pxolist
            Call CopyXOLists(XOSize, TempDone(), TempXOList(), PXOList(), PCurrentXOver(), XOverList(), CurrentXOver(), NumRecsI())
            
        End If
       
    'end section 20*********************************************
    '7.952 5k perms
    '0.360 5k perms
        
        xNextno = Nextno
        'Check for any evidence in the new seqs
            
            Call CheckDrop(Steps(), SEventNumber, StepNo, Nextno, oNextno, NumRecsI(), RedoListSize, RedoList())
            
            X = oNextno + 1
            
           
            Call DropSeqs(X, StepNo, Nextno, Steps(), RedoListSize, RedoList(), StraiName(), MissingData(), SeqNum(), PermValid(), PermDiffs(), TempDone(), PXOList(), PCurrentXOver(), TraceSub(), ActualSeqSize(), MinSeqSize, NumRecsI())

        
        If Nextno > -1 Then
            If Nextno < xNextno Then
                
                
                ReDim TValid(Nextno, Nextno), TDiffs(Nextno, Nextno)
                For X = 0 To Nextno
                    For Y = 0 To Nextno
                        TValid(X, Y) = PermValid(X, Y)
                        TDiffs(X, Y) = PermDiffs(X, Y)
                                
                    Next Y
                Next X
                ReDim PermValid(Nextno, Nextno), PermDiffs(Nextno, Nextno)
                For X = 0 To Nextno
                    For Y = 0 To Nextno
                        PermValid(X, Y) = TValid(X, Y)
                        PermDiffs(X, Y) = TDiffs(X, Y)
                                
                    Next Y
                Next X
            End If
            
            ReDim DoneSeq(Nextno, UBound(PXOList, 2))
            For X = 0 To Nextno
                For Y = 0 To PCurrentXOver(X)
                    If Y <= UBound(TempDone, 2) Then
                        DoneSeq(X, Y) = TempDone(X, Y)
                    End If
                Next Y
            Next X
        End If
        Eventnumber = Eventnumber + EventAdd
    
    If SNextno < Nextno Then
        ReDim Relevant(Nextno), UnInvolved(Nextno)
        
        ReDim TempXOList(Nextno, UBound(PXOList, 2))
        TXOS = Nextno
    End If
    
    
    
    
Loop



End Sub
Public Sub TestMoveInTree(BootFlag As Byte, BPos3 As Long, EPos3 As Long, SeqPair() As Byte, MinPair() As Byte, ISeqs() As Long, SeqNum() As Integer)
Dim ValidX() As Double, DiffsX() As Double
Dim Reps As Integer, tMat() As Double, CutOff As Double, tFCMat() As Double, tSCMat() As Double, TraceBak() As Double, OS As String
Dim TotAdd As Long
Dim DLen() As Double, LLen() As Double, Treebyte() As Byte, NodeDepth() As Integer, BootDepth() As Integer, ReplaceVal As Double, AvDX As Double, K As Long
Dim fMatAv() As Double, fTreeAv() As Double, sMatAv() As Double, sTreeAv() As Double
Dim EraseF As Byte
Dim tSMat() As Double, tFMat() As Double, MinDist(1) As Double, Outlyer(2) As Byte, tSeqNumF() As Integer, tSeqNumS() As Integer
Dim Valtot(1) As Long, MaxMD() As Long, Weight() As Long, Location() As Long, Ally() As Long, Alias() As Long, Px() As Integer, xx1() As Integer, xx2() As Integer
Dim LTree(1) As Long, FHolder() As Byte, SHolder() As Byte, Prod1() As Double, Prod2() As Double, Prod3() As Double, DEN() As Long
Dim UB As Long, NSeqs As Long, TraceSeqs() As Long, tFAMat() As Double, tSAMat() As Double
Dim SCO As Long


CutOff = 50
If DoQuick = 0 Then
    Reps = 10
Else
    Reps = 0
End If
If BPos3 < EPos3 Then
    SCO = CLng((EPos3 - BPos3) / 2)
Else
    SCO = CLng((EPos3 + Len(StrainSeq(0)) - BPos3) / 2)
End If
If SCO > 20 Then SCO = 20

If LongWindedFlag = 1 Then
    mp0 = -1
    If SEventNumber > 0 Then
        If Nextno <> UBound(FMat, 1) Then
            ReDim FMat(Nextno, Nextno)
            ReDim SMat(Nextno, Nextno)
            ReDim SubValid(Nextno, Nextno)
            ReDim SubDiffs(Nextno, Nextno)
        End If
    Else
        ReDim FMat(Nextno, Nextno)
        ReDim SMat(Nextno, Nextno)
        ReDim SubValid(Nextno, Nextno)
        ReDim SubDiffs(Nextno, Nextno)
    End If
    UB = UBound(PermValid, 1)
    efl = 0
    
    If LongFlag = 0 Then
        Dummy = vQuickDist(Len(StrainSeq(0)), Nextno, UB, BPos3, EPos3, FMat(0, 0), SMat(0, 0), SubValid(0, 0), SubDiffs(0, 0), PermValid(0, 0), PermDiffs(0, 0), SeqNum(0, 0), ISeqs(0))
        MinDist(0) = 1000000
        MinDist(1) = 1000000
        
        Outlyer(0) = 2
        Outlyer(1) = 1
        Outlyer(2) = 0
        z = 0
            For X = 0 To 1
                For Y = X + 1 To 2
                    If FMat(ISeqs(X), ISeqs(Y)) < MinDist(0) Then
                        MinDist(0) = FMat(ISeqs(X), ISeqs(Y))
                        MinPair(0) = z
                        SeqPair(0) = X
                        SeqPair(1) = Y
                        SeqPair(2) = Outlyer(z)
                    End If
                    
                    If SMat(ISeqs(X), ISeqs(Y)) < MinDist(1) Then
                        MinDist(1) = SMat(ISeqs(X), ISeqs(Y))
                        MinPair(1) = z
                    End If
                    z = z + 1
                Next Y
            Next X
            
            If MinPair(0) = MinPair(1) Then
                ReDim FCMat(Nextno, Nextno), SCMat(Nextno, Nextno)
                Exit Sub
                efl = 1
            End If
    End If
    
    '0.021
    
    
    Dummy = QuickDist(Len(StrainSeq(0)), Nextno, UB, BPos3, EPos3, FMat(0, 0), SMat(0, 0), SubValid(0, 0), SubDiffs(0, 0), PermValid(0, 0), PermDiffs(0, 0), SeqNum(0, 0))
    
    
    
    
    '0.240
    
   
    
    'WHat is the situation?
    'Which pair is most closely related in the backgound?
    For X = 0 To Nextno
        FMat(X, X) = 0
    Next X
    If X = 12345 Then
        Open "distmat.csv" For Output As #1
        
        For X = 0 To Nextno
            OS = ""
                
            For Y = 0 To Nextno
                If PermValid(X, Y) > 0 Or X = X Then
                    OS = OS + " 0" + Trim(Str(CLng(FMat(X, Y) * 10000) / 10000))
                Else
                    OS = OS + " 0.0000"
                End If
            Next Y
            Print #1, OS
        Next X
        Close #1
        
        X = X
    End If

Else
    ReDim xx1(3), xx2(3), Prod1(Len(StrainSeq(0))), Prod2(Len(StrainSeq(0))), Prod3(Len(StrainSeq(0))), Alias(Len(StrainSeq(0))), Ally(Len(StrainSeq(0))), Location(Len(StrainSeq(0))), Weight(0, Len(StrainSeq(0)))
    ReDim tSeqNumF(Len(StrainSeq(0)), Nextno), tSeqNumS(Len(StrainSeq(0)), Nextno)
    If BPos3 < EPos3 Then
        For X = 1 To BPos3 - 1
            For z = 0 To Nextno
                tSeqNumF(X, z) = SeqNum(X, z)
                tSeqNumS(X, z) = 46
            Next z
        Next X
        For X = BPos3 To EPos3
            For z = 0 To Nextno
                tSeqNumS(X, z) = SeqNum(X, z)
                tSeqNumF(X, z) = 46
            Next z
        Next X
        For X = EPos3 + 1 To Len(StrainSeq(0))
            For z = 0 To Nextno
                tSeqNumF(X, z) = SeqNum(X, z)
                tSeqNumS(X, z) = 46
            Next z
        Next X
    Else
        For X = 1 To EPos3
            For z = 0 To Nextno
                tSeqNumS(X, z) = SeqNum(X, z)
                tSeqNumF(X, z) = 46
            Next z
        Next X
        For X = EPos3 + 1 To BPos3 - 1
            For z = 0 To Nextno
                tSeqNumF(X, z) = SeqNum(X, z)
                tSeqNumS(X, z) = 46
            Next z
        Next X
        For X = BPos3 To Len(StrainSeq(0))
            For z = 0 To Nextno
                tSeqNumS(X, z) = SeqNum(X, z)
                tSeqNumF(X, z) = 46
            Next z
        Next X
    End If
    
    'ee = GetTickCount
    'tt = ee - ss
    'X = X
    
    ReDim FMat(Nextno, Nextno)
    ReDim Px(Nextno, Len(StrainSeq(0)))
    DNADIST 1, 0.5, 0, 0, 0.25, 0.25, 0.25, 0.25, Nextno + 1, Len(StrainSeq(0)), tSeqNumF(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), FMat(0, 0)
    ReDim SMat(Nextno, Nextno)
    ReDim Px(Nextno, Len(StrainSeq(0)))
    '
    DNADIST 1, 0.5, 0, 0, 0.25, 0.25, 0.25, 0.25, Nextno + 1, Len(StrainSeq(0)), tSeqNumS(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), SMat(0, 0)
End If
 

For X = 0 To Nextno
    For Y = X + 1 To Nextno
        FMat(X, Y) = (CLng(FMat(X, Y) * 10000000)) / 10000000
        FMat(Y, X) = FMat(X, Y)
        SMat(X, Y) = (CLng(SMat(X, Y) * 10000000)) / 10000000
        SMat(Y, X) = SMat(X, Y)
    Next Y
Next X
'0.240

'Impute missing entries in distance matrices (ie with pernvalid/subvalid= 0)
'using a nearest neghbours method
If X = 12345 And SEventNumber = 9 Then
        Open "distmat.csv" For Output As #1
        XX = CurDir
        For X = 0 To PermNextNo ' To 0 Step -1  'XXXX0 To
            OS = ""
                
            For Y = 0 To PermNextNo    'xxxx 0 To
                If SubValid(X, Y) > 0 Or X = X Then
                    OS = OS + " 0" + Trim(Str(CLng(FMat(X, Y) * 10000) / 10000))
                Else
                    OS = OS + " 0.0000"
                End If
            Next Y
            Print #1, OS
        Next X
        Close #1
        X = X
    End If
If LongWindedFlag = 1 Then

    'Dummy = CheckMatrix(UB, SCO, MinSeqSize, Nextno, ValTot(0), ISeqs(0), PermValid(0, 0), SubValid(0, 0), FMat(0, 0), SMat(0, 0), FAMat(0, 0), SAMat(0, 0))
    Call CheckMatrixX(UB, SCO, MinSeqSize, Nextno, Valtot(), ISeqs(), PermValid(), SubValid(), FMat(), SMat(), FAMat(), SAMat())

End If

 '0.060

'Make bAckup of Fmat (it is destroyed by neighbour)
ReDim TraceSeqs(1, Nextno)

NSeqs = -1
For X = 0 To Nextno
    If FMat(X, X) <> 3 Then
        NSeqs = NSeqs + 1
        TraceSeqs(0, X) = NSeqs
        TraceSeqs(1, NSeqs) = X
    
    End If
Next X
'0.260


If NSeqs > 2 Then
    ReDim tFMat(NSeqs, NSeqs), tSMat(NSeqs, NSeqs), tFAMat(NSeqs, NSeqs), tSAMat(NSeqs, NSeqs), SHolder((NSeqs + 1) * 40 * 2), FHolder((NSeqs + 1) * 40 * 2)
Else
    MinPair(0) = MinPair(1)
    ReDim FCMat(Nextno, Nextno), SCMat(Nextno, Nextno)
    Exit Sub
End If
A = 0
B = 0
'0.260

For X = 0 To Nextno
    
    If FMat(X, X) <> 3 Then
        tFMat(A, A) = 0
        tSMat(A, A) = 0
        B = A + 1
        
        For Y = X + 1 To Nextno
            If FMat(Y, Y) <> 3 Then
                tFMat(A, B) = FMat(X, Y)
                tSMat(A, B) = SMat(X, Y)
                tFMat(B, A) = FMat(X, Y)
                tSMat(B, A) = SMat(X, Y)
                B = B + 1
            End If
            
        Next Y
        A = A + 1
   
    End If
Next X

'0.270
'find outie in fmat
    Dim TotMat() As Double, MaxSc As Double
    ReDim TotMat(NSeqs)
    For X = 0 To NSeqs
        For Y = X + 1 To NSeqs
            If tFMat(X, X) < 3 And tFMat(Y, Y) < 3 Then
                TotMat(X) = TotMat(X) + tFMat(X, Y)
                TotMat(Y) = TotMat(Y) + tFMat(X, Y)
            End If
        Next Y
    Next X
    MaxSc = 0
    For X = 0 To NSeqs
        If TotMat(X) > MaxSc Then
            Outie = X
            MaxSc = TotMat(X)
        End If
    Next X
'0.260

If X = 12345 And SEventNumber = 66 Then
        Open "distmat.csv" For Output As #1
        zzzzzz = CurDir
        For X = 0 To NSeqs ' To 0 Step -1  'XXXX0 To
            OS = ""
                
            For Y = 0 To NSeqs    'xxxx 0 To
                OS = OS + Str(tFMat(X, Y)) + ","
                
            Next Y
            Print #1, OS
        Next X
        Close #1
        X = X
    End If

ReDim ColTotals(NSeqs)
'XX = UBound(FHolder, 1) '1,0,61,64,63,5120,63
LTree(0) = NEIGHBOUR(1, 0, TRndSeed, Outie + 1, NSeqs + 1, tFMat(0, 0), FHolder(0), ColTotals(0), tFAMat(0, 0))
ReDim ColTotals(NSeqs)
LTree(1) = NEIGHBOUR(1, 0, TRndSeed, Outie + 1, NSeqs + 1, tSMat(0, 0), SHolder(0), ColTotals(0), tSAMat(0, 0))

'0.290


'Outputs trees if you want these for debugging purposes
If X = 12345 Then
    xxx = CurDir
    ssy = ""
    XX = BPos3
    XX = EPos3
    SX = ""
    
    For z = 1 To LTree(1) '1264
        ssy = ssy + Chr(SHolder(z))
    Next z '

    For z = 1 To LTree(0) '1264
        SX = SX + Chr(FHolder(z))
    Next z '

    Open "testfa" & Str(SEventNumber + 1) & ".tre" For Output As #1
    Print #1, SX
   
    Close #1
     Open "testfb" & Str(SEventNumber + 1) & ".tre" For Output As #1
   
    Print #1, ssy
    Close #1
    X = X
    XX = TraceSeqs(1, 60)
End If

NameLen = Len(Trim$(CStr(NSeqs)))
If NameLen < 2 Then NameLen = 2




Call Tree2Array(NameLen, NSeqs, LTree(0), FHolder(), tFAMat())
Call Tree2Array(NameLen, NSeqs, LTree(1), SHolder(), tSAMat())
'0.311


ReDim FAMat(Nextno, Nextno), SAMat(Nextno, Nextno)
For X = 0 To NSeqs
    For Y = X + 1 To NSeqs
        tFAMat(X, Y) = CLng(tFAMat(X, Y) * 10000) / 10000
        tSAMat(X, Y) = CLng(tSAMat(X, Y) * 10000) / 10000
        tFAMat(Y, X) = tFAMat(X, Y)
        tSAMat(Y, X) = tSAMat(X, Y)
        FAMat(TraceSeqs(1, X), TraceSeqs(1, Y)) = tFAMat(X, Y)
        FAMat(TraceSeqs(1, Y), TraceSeqs(1, X)) = tFAMat(X, Y)
        SAMat(TraceSeqs(1, X), TraceSeqs(1, Y)) = tSAMat(X, Y)
        SAMat(TraceSeqs(1, Y), TraceSeqs(1, X)) = tSAMat(X, Y)
    Next Y
Next X
'0.320

For X = 0 To Nextno
    If FMat(X, X) = 3 Then
        For Y = 0 To Nextno
            FAMat(X, Y) = ((Nextno * 3) - 1) / 1000
            FAMat(Y, X) = ((Nextno * 3) - 1) / 1000
            SAMat(X, Y) = ((Nextno * 3) - 1) / 1000
            SAMat(Y, X) = ((Nextno * 3) - 1) / 1000
        Next Y
    End If
Next X
'0.330




If LongWindedFlag = 0 Or X = X Then
    MinDist(0) = 1000000
    MinDist(1) = 1000000
    
    Outlyer(0) = 2
    Outlyer(1) = 1
    Outlyer(2) = 0
    z = 0
    
    
    'Eventnumber = Eventnumber
    'Nextno = Nextno
    For X = 0 To 1
        For Y = X + 1 To 2
            If FAMat(ISeqs(X), ISeqs(Y)) < MinDist(0) Then '15-4 = 0.045: 15-11 = 0.042, 4-11 = 0.045
                MinDist(0) = FAMat(ISeqs(X), ISeqs(Y))
                MinPair(0) = z
                SeqPair(0) = X
                SeqPair(1) = Y
                SeqPair(2) = Outlyer(z)
            End If
            If SAMat(ISeqs(X), ISeqs(Y)) < MinDist(1) Then '15-4 = 0.04, 15-11 = 0.014, 4-11=0.04
                MinDist(1) = SAMat(ISeqs(X), ISeqs(Y))
                MinPair(1) = z
            End If
            z = z + 1
        Next Y
    Next X
End If

'0.320
If Reps > 0 And BootFlag = 1 And (MinPair(1) <> MinPair(0) Or LongFlag = 1) Then
    
    Dim TMatch() As Byte
    
    Form1.SSPanel1.Caption = "Making bootstrapped NJ tree"
    ODir = CurDir
    ChDir App.Path
    ChDrive App.Path
    'Perform bootstrap replicates
    Dim tSeqNum2() As Integer, tSeqNum() As Integer, WeightMod() As Long, Scratch() As Integer, Length As Long
    Dim FF As Long, TreeString As String, DstMat() As Double, Num1() As Long, Num2() As Long, Num() As Double
    ReDim tSeqNum2(Len(StrainSeq(0)), NSeqs)
    
    
    Dummy = MaketSeqNum(Len(StrainSeq(0)), Nextno, tSeqNum2(0, 0), SeqNum(0, 0), FMat(0, 0))
    
    
    
    'Call MakeConsenseFiles(NSeqs)
    '0.411
    
    For X = 0 To 1
        FF = FreeFile
        
        'Close #FF
        
        If X = 0 Then 'ie inner alignment
            Call MakeETSeqNum(NSeqs, Length, BPos3, EPos3, tSeqNum(), tSeqNum2())
            Call MakeNodeDepth(NSeqs, TraceBak(), tSAMat(), NodeDepth())
        Else
            Call MakeETSeqNum(NSeqs, Length, EPos3 + 1, BPos3 - 1, tSeqNum(), tSeqNum2())
            Call MakeNodeDepth(NSeqs, TraceBak(), tFAMat(), NodeDepth())
        End If
        
        
        
        Call MakeTMatch(NSeqs, TMatch(), NodeDepth())
        
        ReDim WeightMod(Reps, Length - 1), Scratch(Length)
        
        '0.461
        '0.441
        
        
        'Dummy = SEQBOOT(BSRndNumSeed, Reps, Length, Scratch(0), WeightMod(0, 0))
        
        Dummy = SEQBOOT2(BSRndNumSeed, Reps, Length, Scratch(0), WeightMod(0, 0))
        
        '0.671
        '0.641
        
        
        'ReDim Alias(Length), Ally(Length), Num(Reps + 1), DEN(Reps + 1), Num1(Reps + 1), Num2(Reps + 1), xx1(3), xx2(3), Weight(Reps, Length), Location(Length), Px(NextNo, Length), Prod1(Length), Prod2(Length), Prod3(Length) 'doub
        
        '0.671
        '0.641
        'Dummy = BootDist(Reps, TCoeffVar, 0.5, 0, NSeqs + 1, Length, tSeqNum(0, 0), Alias(0), Ally(0), Weight(0, 0), Location(0), Px(0, 0), xx1(0), xx2(0), Prod1(0), Prod2(0), Prod3(0), DstMat(0, 0, 0), DistVal(0), Num1(0), Num2(0), DEN(0), Num(0), WeightMod(0, 0))

        ReDim ValidX(Reps), DiffsX(Reps), DstMat(Reps, NSeqs, NSeqs)
        Dummy = FastBootDist(1, Reps, NSeqs, Length, DiffsX(0), ValidX(0), WeightMod(0, 0), tSeqNum(0, 0), DstMat(0, 0, 0))
        'Exit Sub
        '2.864
        '2.153 with fastbootdist
        '2.144,1.933, 1.913, 1.883, 1.882
        If AbortFlag = 1 Then Exit Sub
        
        
        
        
        
        
        ReDim tFMat(NSeqs, NSeqs), tMat(NSeqs, NSeqs), DLen(NSeqs)
        For Y = 1 To Reps
            
            ReDim FHolder(Nextno * 40 * 2)
            Dummy = TransferDist(NSeqs, Y, Reps, tFMat(0, 0), DstMat(0, 0, 0))
            'ReDim ColTotals(NSeqs)
            LTree(0) = NEIGHBOUR(1, 0, TRndSeed, Outie + 1, NSeqs + 1, tFMat(0, 0), FHolder(0), ColTotals(0), tFAMat(0, 0))
            '2.965
           
            Call TreeGroups(NSeqs, FHolder(), LTree(0), NameLen, TMatch(), DLen())
            '2.984
             'Exit Sub
            
            If AbortFlag = 1 Then
                Close #FF
                Exit Sub
            End If
        Next Y
        
        '3.555
        '3.075
        
        For z = 0 To NSeqs
             DLen(z) = CLng((DLen(z) + 1) / (Reps + 1) * 100)
             'X = X '100,100,100,100,82,82,100,100,64,100,100,100,55
        Next z
        'Collapse nodes in FAMat/SAMat with no support - Put these in FCMat and SCMat
        'famat/samat = tree distance matrices - not path lengths but an encoding of
        'the rooted tree topology in a distance matrix
        'fcmat/scmat = tree distance matrices with nodes collapsed
        'Dlen = array containing bootstrap support for nodes
        'tracebak = contains an encoding of the tree topology - ie tells you which node corrsponds to which
        'distance in samat/famat
        If X = 0 Then
            Call CollapseNodes(NSeqs, CutOff, DLen(), TraceBak(), tSAMat(), tSCMat())
        Else
            Call CollapseNodes(NSeqs, CutOff, DLen(), TraceBak(), tFAMat(), tFCMat())
        End If
        '3.075
       
        
        
    Next X
    'Exit Sub
    '9.794
    '8.813
     'Exit Sub
    
    ReDim FCMat(Nextno, Nextno), SCMat(Nextno, Nextno)
    If AbortFlag = 1 Then
        Exit Sub
    End If
    For X = 0 To NSeqs
        For Y = X + 1 To NSeqs
            SCMat(TraceSeqs(1, X), TraceSeqs(1, Y)) = tSCMat(X, Y)
            SCMat(TraceSeqs(1, Y), TraceSeqs(1, X)) = tSCMat(X, Y)
            FCMat(TraceSeqs(1, X), TraceSeqs(1, Y)) = tFCMat(X, Y)
            FCMat(TraceSeqs(1, Y), TraceSeqs(1, X)) = tFCMat(X, Y)
            
        Next Y
    Next X
    
    For X = 0 To Nextno
        If FMat(X, X) = 3 Then
            For Y = 0 To Nextno
                FCMat(X, Y) = ((Nextno * 3) - 1) / 1000
                FCMat(Y, X) = FCMat(X, Y)
                SCMat(X, Y) = ((Nextno * 3) - 1) / 1000
                SCMat(Y, X) = SCMat(X, Y)
            Next Y
        End If
        
    Next X
    For X = 0 To Nextno
        For Y = X + 1 To Nextno
            FCMat(X, Y) = (CLng(FCMat(X, Y) * 10000000)) / 10000000
            FCMat(Y, X) = FCMat(X, Y)
            SCMat(X, Y) = (CLng(SCMat(X, Y) * 10000000)) / 10000000
            SCMat(Y, X) = SCMat(X, Y)
        Next Y
    Next X
    On Error Resume Next
    ChDir ODir
    ChDrive ODir
    On Error GoTo 0
Else
    ReDim FCMat(Nextno, Nextno)
    ReDim SCMat(Nextno, Nextno)
    For X = 0 To Nextno
        For Y = 0 To Nextno
            FCMat(X, Y) = FAMat(X, Y)
            SCMat(X, Y) = SAMat(X, Y)
        Next Y
    Next X
End If


'Make sure diagonals = 0
For X = 0 To Nextno
    FMat(X, X) = 0
    SMat(X, X) = 0
    FAMat(X, X) = 0
    SAMat(X, X) = 0
Next X


End Sub
Public Sub MakeTMatch(Nextno As Long, TMatch() As Byte, NodeDepth() As Integer)
Dim X As Long
ReDim TMatch(Nextno, Nextno)
For X = 0 To Nextno
    For Y = 0 To Nextno
        If NodeDepth(X, Y) > -1 Then
            TMatch(X, NodeDepth(X, Y)) = 1
        Else
            Exit For
        End If
    Next Y
Next X
End Sub

Public Sub TreeGroups(Nextno As Long, THolder() As Byte, TLen As Long, NLen, TMatch() As Byte, DLen() As Double)
Dim X As Long, Y As Long, z As Long, Cnt As Long, SeqID As Long, TArray() As Byte, NCount As Long
ReDim TArray(Nextno, Nextno)
NCount = -1
Dim Miss As Long, Hit As Long, DoneNode() As Byte
ReDim DoneNode(Nextno)

If X = X Then
    Dummy = TreeGroupsX(Nextno, THolder(0), TLen, NLen, DoneNode(0), TArray(0, 0), TMatch(0, 0), DLen(0))
Else
    For X = 1 To TLen
        If THolder(X) = 40 Then 'ie (
            Cnt = 1
            NCount = NCount + 1
            Y = X + 1
            Do While Cnt > 0
                If THolder(Y) = 40 Then 'ie (
                    Cnt = Cnt + 1
                ElseIf THolder(Y) = 41 Then 'ie )
                    Cnt = Cnt - 1
                ElseIf THolder(Y) = 83 Then 'ie S
                    SeqID = 0
                    For z = 1 To NLen
                        SeqID = SeqID + CLng(THolder(Y + z) - 48) * 10 ^ (NLen - z)
                    Next z
                    TArray(NCount, SeqID) = 1
                End If
                Y = Y + 1
            Loop
        End If
    
    Next X
    
    
    
    For X = 0 To Nextno
        For Y = 0 To Nextno
            If DoneNode(Y) = 0 Then
                Miss = 0
                Hit = 0
                For z = 0 To Nextno
                    If TArray(Y, z) = TMatch(X, z) Then
                        Hit = Hit + 1
                    Else
                        Miss = Miss + 1
                    End If
                    If Miss > 0 And Hit > 0 Then Exit For
                Next z
                If Miss = 0 Or Hit = 0 Then
                    DLen(X) = DLen(X) + 1
                    DoneNode(Y) = 1
                Else
                    X = X
                End If
            End If
        Next Y
    Next X
    X = X
End If
End Sub
