Attribute VB_Name = "CmdOutput"
'option explicit
''''''''''''''''''''''''''''''''''''''''
' Joacim Andersson, Brixoft Software
' http://www.brixoft.net
''''''''''''''''''''''''''''''''''''''''

' STARTUPINFO flags
Private Const STARTF_USESHOWWINDOW = &H1
Private Const STARTF_USESTDHANDLES = &H100

' ShowWindow flags
Private Const SW_HIDE = 0

' DuplicateHandle flags
Private Const DUPLICATE_CLOSE_SOURCE = &H1
Private Const DUPLICATE_SAME_ACCESS = &H2

' Error codes
Private Const ERROR_BROKEN_PIPE = 109

Private Type SECURITY_ATTRIBUTES
    nLength As Long
    lpSecurityDescriptor As Long
    bInheritHandle As Long
End Type

Private Type STARTUPINFO
    CB As Long
    lpReserved As String
    lpDesktop As String
    lpTitle As String
    dwX As Long
    dwY As Long
    dwXSize As Long
    dwYSize As Long
    dwXCountChars As Long
    dwYCountChars As Long
    dwFillAttribute As Long
    dwFlags As Long
    wShowWindow As Integer
    cbReserved2 As Integer
    lpReserved2 As Long
    hStdInput As Long
    hStdOutput As Long
    hStdError As Long
End Type

Private Type PROCESS_INFORMATION
    hProcess As Long
    hThread As Long
    dwProcessId As Long
    dwThreadId As Long
End Type

Private Declare Function CreatePipe _
 Lib "kernel32" ( _
 phReadPipe As Long, _
 phWritePipe As Long, _
 lpPipeAttributes As Any, _
 ByVal nSize As Long) As Long

Private Declare Function ReadFile _
 Lib "kernel32" ( _
 ByVal hFile As Long, _
 lpBuffer As Any, _
 ByVal nNumberOfBytesToRead As Long, _
 lpNumberOfBytesRead As Long, _
 lpOverlapped As Any) As Long

Private Declare Function CreateProcess _
 Lib "kernel32" Alias "CreateProcessA" ( _
 ByVal lpApplicationName As String, _
 ByVal lpCommandLine As String, _
 lpProcessAttributes As Any, _
 lpThreadAttributes As Any, _
 ByVal bInheritHandles As Long, _
 ByVal dwCreationFlags As Long, _
 lpEnvironment As Any, _
 ByVal lpCurrentDriectory As String, _
 lpStartupInfo As STARTUPINFO, _
 lpProcessInformation As PROCESS_INFORMATION) As Long

Private Declare Function GetCurrentProcess _
 Lib "kernel32" () As Long

Private Declare Function DuplicateHandle _
 Lib "kernel32" ( _
 ByVal hSourceProcessHandle As Long, _
 ByVal hSourceHandle As Long, _
 ByVal hTargetProcessHandle As Long, _
 lpTargetHandle As Long, _
 ByVal dwDesiredAccess As Long, _
 ByVal bInheritHandle As Long, _
 ByVal dwOptions As Long) As Long

Private Declare Function CloseHandle _
 Lib "kernel32" ( _
 ByVal hObject As Long) As Long

Private Declare Function OemToCharBuff _
 Lib "user32" Alias "OemToCharBuffA" ( _
 lpszSrc As Any, _
 ByVal lpszDst As String, _
 ByVal cchDstLength As Long) As Long

' Function GetCommandOutput
'
' sCommandLine:  [in] Command line to launch
' fStdOut        [in,opt] True (defualt) to capture output to STDOUT
' fStdErr        [in,opt] True to capture output to STDERR. False is default.
' fOEMConvert:   [in,opt] True (default) to convert DOS characters to Windows, False to skip conversion
'
' Returns:       String with STDOUT and/or STDERR output
'
Public Function GetCommandOutput(sCommandLine As String, PrF, Optional fStdOut As Boolean = True, _
                                 Optional fStdErr As Boolean = False, Optional fOEMConvert As Boolean = True) As String

    Dim Pos As Long, LastPos As Long
    Dim hPipeRead As Long, hPipeWrite1 As Long, hPipeWrite2 As Long, SS As Long, EE As Long, TT As Long, LV As Double
    Dim hCurProcess As Long
    Dim sa As SECURITY_ATTRIBUTES
    Dim si As STARTUPINFO
    Dim pi As PROCESS_INFORMATION
    Dim baOutput() As Byte
    Dim sNewOutput As String
    Dim lBytesRead As Long, SCount As Long, FLen As Long, LastLen As Long, CurPos As Long, NewPos As Long, Count As Long, NewSurface As Long
    Dim fTwoHandles As Boolean
    Dim Pict As Long
    Dim lRet As Long
    Dim PntAPI As POINTAPI
    Dim x As Long
    Dim Pos2 As Long
    Dim TS As String
    Dim BUFSIZE As Long
    DoingShellFlag = DoingShellFlag + 1
    If PrF = 2 Then
        BUFSIZE = 1024
    Else
        BUFSIZE = 1024      ' pipe buffer size
    End If
    ' At least one of them should be True, otherwise there's no point in calling the function
    If (Not fStdOut) And (Not fStdErr) Then
        Err.Raise 5         ' Invalid Procedure call or Argument
    End If

    ' If both are true, we need two write handles. If not, one is enough.
    fTwoHandles = fStdOut And fStdErr

    ReDim baOutput(BUFSIZE - 1) As Byte

    With sa
        .nLength = Len(sa)
        .bInheritHandle = 1    ' get inheritable pipe handles
    End With

    If CreatePipe(hPipeRead, hPipeWrite1, sa, BUFSIZE) = 0 Then
        DoingShellFlag = DoingShellFlag - 1
        Exit Function
    End If

    hCurProcess = GetCurrentProcess()

    ' Replace our inheritable read handle with an non-inheritable. Not that it
    ' seems to be necessary in this case, but the docs say we should.
    Call DuplicateHandle(hCurProcess, hPipeRead, hCurProcess, hPipeRead, 0&, _
                         0&, DUPLICATE_SAME_ACCESS Or DUPLICATE_CLOSE_SOURCE)

    ' If both STDOUT and STDERR should be redirected, get an extra handle.
    If fTwoHandles Then
        Call DuplicateHandle(hCurProcess, hPipeWrite1, hCurProcess, hPipeWrite2, 0&, _
                             1&, DUPLICATE_SAME_ACCESS)
    End If

    With si
        .CB = Len(si)
        .dwFlags = STARTF_USESHOWWINDOW Or STARTF_USESTDHANDLES
        .wShowWindow = SW_HIDE          ' hide the window

        If fTwoHandles Then
            .hStdOutput = hPipeWrite1
            .hStdError = hPipeWrite2
        ElseIf fStdOut Then
            .hStdOutput = hPipeWrite1
        Else
            .hStdError = hPipeWrite1
        End If
    End With

    If CreateProcess(vbNullString, sCommandLine, ByVal 0&, ByVal 0&, 1, 0&, _
     ByVal 0&, vbNullString, si, pi) Then

        ' Close thread handle - we don't need it
        Call CloseHandle(pi.hThread)

        ' Also close our handle(s) to the write end of the pipe. This is important, since
        ' ReadFile will *not* return until all write handles are closed or the buffer is full.
        Call CloseHandle(hPipeWrite1)
        hPipeWrite1 = 0
        If hPipeWrite2 Then
            Call CloseHandle(hPipeWrite2)
            hPipeWrite2 = 0
        End If
        Dim LTM As Variant, oLTM As Variant, FFX As Long, LFX As Long
        LTM = 0
        Do
            ' Add a DoEvents to allow more data to be written to the buffer for each call.
            ' This results in fewer, larger chunks to be read.
            
            'DoEvents
            
            'If PrF <> 2 Then
                If LTM < 10000 Or PrF <> 14 Then 'checks for a pipe blockage but will only do something different with phml
                    oLTM = LTM
                    LTM = Abs(GetTickCount)
                    If ReadFile(hPipeRead, baOutput(0), BUFSIZE, lBytesRead, ByVal 0&) = 0 Then
                        Exit Do
                    End If
                    If AbortFlag = 1 Then
                        Exit Do
                    End If
                    LTM = Abs(GetTickCount) - LTM
                    If LTM < 10 Then
                        LTM = oLTM * 2
                    End If
                    
                
                ElseIf PrF = 14 Then
                    lBytesRead = 0
                    'baOutput(0) = 0
                    FFX = FreeFile
                    LFX = 0
                    If DebuggingFlag < 2 Then On Error Resume Next
                    Open "treefile" For Input As FFX
                    LFX = LOF(FFX)
                    Close FFX
                    On Error GoTo 0
                    If LFX > 0 Then
                        Exit Do
                    End If
                    LTM = LTM - 100
                End If
                
                DoEvents
            'Else
            '    lBytesRead = 0
            'End If
            'If lBytesRead = 0 > 0 Then
                If fOEMConvert Then
                    ' convert from "DOS" to "Windows" characters
                    sNewOutput = String$(lBytesRead, 0)
                    Call OemToCharBuff(baOutput(0), sNewOutput, lBytesRead)
                    'X = X
                Else
                    ' perform no conversion (except to Unicode)
                    sNewOutput = Left$(StrConv(baOutput(), vbUnicode), lBytesRead)
                End If
            If DebuggingFlag < 2 Then On Error Resume Next
    
                GetCommandOutput = GetCommandOutput & sNewOutput
            'End If
                LV = Len(GetCommandOutput)
            On Error GoTo 0
            'Open "screen.out" For Output As #1
            'Print #1, GetCommandOutput
            'Close #1
            ' If you are executing an application that outputs data during a long time,
            ' and don't want to lock up your application, it might be a better idea to
            ' wrap this code in a class module in an ActiveX EXE and execute it asynchronously.
            ' Then you can raise an event here each time more data is available.
            'RaiseEvent OutputAvailabele(sNewOutput)
            'Form1.SSPanel1.Caption = Right$(GetCommandOutput, 50)
            SS = Abs(GetTickCount)
            DoEvents
            If Abs(SS - EE) > 500 Then
                            DoEvents
                EE = SS
                If PrF = 23 Then ' interval
                
                    TS = Right$(GetCommandOutput, 200)
                    LastPos = InStr(1, TS, "Run", vbBinaryCompare)
                    If LastPos > 0 Then
                        Pos = InStr(LastPos + 1, TS, ":", vbBinaryCompare)
                        If Pos > 0 Then
                            TS = Trim(Mid$(TS, LastPos + 3, Pos - LastPos - 3))
                            If DebuggingFlag < 2 Then On Error Resume Next
                            TS = Trim(Str(CLng(val(TS))))
                            Form1.SSPanel1.Caption = Trim(Str(TS)) + " of " & Trim(Str(FullSize)) & "  MCMC updates completed"
                            Form1.SSPanel1.Refresh
                            On Error GoTo 0
                        End If
                    End If
                    x = x
                ElseIf PrF = 1 Or PrF = 2 Then 'lard with 1 breakpoint
                    
                    Count = 0
                    Pos = InStr(1, GetCommandOutput, "[.", vbBinaryCompare)
                    
                    Do While Pos > 0
                        LastPos = Pos + 1
                        Pos = InStr(LastPos, GetCommandOutput, ".", vbBinaryCompare)
                        Count = Count + 1
                        If Pos - LastPos > 5 Then Exit Do
                    Loop
                    LV = Form1.ProgressBar1.Value
                    If Count > 0 Then
                        If Count <= 22 Then
                            Form1.ProgressBar1.Value = 5 + (Count / 22) * 90
                        Else
                            Form1.ProgressBar1.Value = 90
                        End If
                        Call UpdateF2Prog
                    End If
                    
                    
                    If LV <> Form1.ProgressBar1.Value Then
                        Form1.Refresh
                    End If
                
                ElseIf PrF = 3 Then
                    LV = Len(GetCommandOutput) / ExpectFL
                    Open "screen.out" For Output As #1
                    Print #1, GetCommandOutput
                    Close #1
                    If LV < 1 Then
                        Form1.ProgressBar1 = LV * 100
                        Form1.SSPanel1.Caption = CStr(CLng(Len(GetCommandOutput) / 38)) + " of " + CStr(CLng(ExpectFL / 38)) + " windows examined"
                        Call UpdateF2Prog
                    End If
                ElseIf PrF = 5 Then
                    
                    
                    'Open "screen.out" For Output As #1
                    'Print #1, GetCommandOutput
                    'Close #1
                    x = val(Right(GetCommandOutput, 3))
                    'X = Len(GetCommandOutput)
                    If x > 0 Then
                        x = x - ExpectFL
                        Form1.ProgressBar1 = x / (x - ExpectFL)
                        Call UpdateF2Prog
                    Else
                        GetCommandOutput = ""
                    End If
                ElseIf PrF = 7 Then
                    GetCommandOutput = Right(GetCommandOutput, 70)
                    Pos = Len(GetCommandOutput)
                    Pos = InStr(1, GetCommandOutput, "LogLk =", vbBinaryCompare)
                    'TS = CurDir
                    'Open "xx.txt" For Append As #10
                    'Print #10, GetCommandOutput
                    'Print #10, ""
                    'Close #10
                    If Pos > 0 Then
                        'Pos2 = InStr(Pos, GetCommandOutput, " ", vbBinaryCompare)
                        TS = Mid(GetCommandOutput, Pos, 18)
                        Do
                            Pos = InStr(1, TS, "  ", vbBinaryCompare)
                            If Pos > 0 Then
                                TS = Left(TS, Pos) + Right(TS, Len(TS) - (Pos + 1))
                                
                            Else
                                Exit Do
                            End If
                        Loop
                        Form1.SSPanel1.Caption = "Making ML Tree (with FastTree; " + Trim(TS) + ")"
                        Form1.ProgressBar1.Value = 70
                        Call UpdateF2Prog
                    Else
                        Pos = InStr(1, GetCommandOutput, "seconds:", vbBinaryCompare)
                        If Pos > 0 Then
                            
                            
                            
                            If Len(GetCommandOutput) > Pos + 9 Then
                                TS = Mid(GetCommandOutput, Pos + 9, Len(GetCommandOutput) - (Pos + 10))
                                
                                Do
                                    Pos = InStr(1, TS, "  ", vbBinaryCompare)
                                    If Pos > 0 Then
                                        TS = Left(TS, Pos) + Right(TS, Len(TS) - (Pos + 1))
                                        x = x
                                    Else
                                        Exit Do
                                    End If
                                Loop
                                
                                Form1.SSPanel1.Caption = "Making ML Tree (with FastTree; " + Trim(TS) + ")"
                                Call UpdateF2Prog
                            End If
                        End If
                        If Form1.ProgressBar1 < 40 Then
                            Pos = InStr(1, GetCommandOutput, "Optimizing GTR", vbBinaryCompare)
                            If Pos > 0 Then
                                Form1.ProgressBar1.Value = 40
                                Call UpdateF2Prog
                            End If
                        End If
                    End If
                    Sleep 200
                ElseIf PrF = 8 Then
                    GetCommandOutput = Right(GetCommandOutput, 70)
                    Pos = Len(GetCommandOutput)
                    Pos = InStr(1, GetCommandOutput, "LogLk =", vbBinaryCompare)
                    'TS = CurDir
                    'Open "xx.txt" For Append As #10
                    'Print #10, GetCommandOutput
                    'Print #10, ""
                    'Close #10
                    If Pos > 0 Then
                        'Pos2 = InStr(Pos, GetCommandOutput, " ", vbBinaryCompare)
                        TS = Mid(GetCommandOutput, Pos, 18)
                        Do
                            Pos = InStr(1, TS, "  ", vbBinaryCompare)
                            If Pos > 0 Then
                                TS = Left(TS, Pos) + Right(TS, Len(TS) - (Pos + 1))
                                
                            Else
                                Exit Do
                            End If
                        Loop
                        Form1.SSPanel1.Caption = "Making FastNJ Tree (with FastTree; " + Trim(TS) + ")"
                        Call UpdateF2Prog
                    Else
                        Pos = InStr(1, GetCommandOutput, "seconds:", vbBinaryCompare)
                        If Pos > 0 Then
                            
                            
                            
                            If Len(GetCommandOutput) > Pos + 9 Then
                                TS = Mid(GetCommandOutput, Pos + 9, Len(GetCommandOutput) - (Pos + 10))
                                
                                Do
                                    Pos = InStr(1, TS, "  ", vbBinaryCompare)
                                    If Pos > 0 Then
                                        TS = Left(TS, Pos) + Right(TS, Len(TS) - (Pos + 1))
                                        x = x
                                    Else
                                        Exit Do
                                    End If
                                Loop
                                
                                Form1.SSPanel1.Caption = "Making FastNJ Tree (with FastTree; " + Trim(TS) + ")"
                                Call UpdateF2Prog
                            End If
                        End If
                        
                    End If
                    Sleep 200
                ElseIf PrF = 14 Then
                    GetCommandOutput = Right(GetCommandOutput, 10000)
                    Pos = InStr(1, GetCommandOutput, "bootstrap analysis", vbBinaryCompare)
                    'XX = Len(GetCommandOutput)
                    'XX = Right$(GetCommandOutput, 120)
                    If Pos = 0 Then
                        If Len(GetCommandOutput) > 50 Then
                            Pos = InStr(Len(GetCommandOutput) - 50, GetCommandOutput, " -", vbBinaryCompare)
                        End If
                        
                        If Pos > 0 Then
                            Pos2 = InStr(Pos, GetCommandOutput, "]", vbBinaryCompare)
                            If Pos2 > 0 Then
                                TS = Mid(GetCommandOutput, Pos, Pos2 - Pos)
                                If TBSReps > 1 And BSupTest = 0 Then
                                    
                                    If Form1.ProgressBar1.Value > 50 Then Form1.ProgressBar1.Value = 40
                                    Pos = 50 - Form1.ProgressBar1.Value
                                    Form1.ProgressBar1.Value = Form1.ProgressBar1.Value + Pos / 20
                                Else
                                    Pos = 95 - Form1.ProgressBar1.Value
                                    Form1.ProgressBar1.Value = Form1.ProgressBar1.Value + Pos / 20
                                End If
                                
                                Form1.SSPanel1.Caption = "Making ML Tree (with PhyML3; LogLk = " + Trim(TS) + ")"
                                Call UpdateF2Prog
                            End If
                        End If
                    Else
                        Pos = InStr(1, GetCommandOutput, "bootstrap analysis", vbBinaryCompare)
                        If Pos = 0 And Len(GetCommandOutput) > 0 And Left(GetCommandOutput, 1) = "." Then
                            TS = Trim(Str(Len(GetCommandOutput)))
                            Form1.ProgressBar1.Value = 50 + (Pos / TBSReps) * 45
                            TS = TS + "% of bootstrap replicates completed"
                            Form1.SSPanel1.Caption = "Making ML Tree (with PhyML3; " + Trim(TS) + ")"
                            Call UpdateF2Prog
                        ElseIf Pos > 0 Then
                            Pos2 = InStr(1, GetCommandOutput, "[.", vbBinaryCompare)
                            'XX = Right$(GetCommandOutput, 20)
                            If Pos2 > 0 Then
                                'count the dots
                                Pos = 1
                                Do
                                    Pos2 = InStr(Pos2 + 1, GetCommandOutput, ".", vbBinaryCompare)
                                    If Pos2 > 0 Then
                                        Pos = Pos + 1
                                    Else
                                        Exit Do
                                    End If
                                    
                                    
                                Loop
                                'Pos = Len(GetCommandOutput) - Pos2
                                TS = Trim(Str((Pos / TBSReps) * 100))
                                'TS = Trim(Str(Pos))
                                Form1.ProgressBar1.Value = 50 + (Pos / TBSReps) * 45
                                TS = TS + "% of bootstrap replicates completed"
                                Form1.SSPanel1.Caption = "Making ML Tree (with PhyML3; " + Trim(TS) + ")"
                                Call UpdateF2Prog
                            End If
                        Else
                            Pos = InStr(1, GetCommandOutput, "aLRT", vbBinaryCompare)
                            If Pos > 0 Then
                                Form1.SSPanel1.Caption = "Making ML Tree (with PhyML3; Computing aLRT branch supports)"
                                Form1.ProgressBar1.Value = 95
                                Call UpdateF2Prog
                            End If
                        End If
                    End If
                    
                    Sleep 100
                ElseIf PrF = 6 Then
                    GetCommandOutput = Right(GetCommandOutput, 50)
                    Do
                        Pos = InStr(1, GetCommandOutput, Chr(13), vbBinaryCompare)
                        If Pos = 0 Then
                            Pos = InStr(1, GetCommandOutput, Chr(10), vbBinaryCompare)
                            If Pos = 0 Then
                                Exit Do
                            Else
                                Mid$(GetCommandOutput, Pos, 1) = " "
                            End If
                        Else
                            x = x
                            Mid$(GetCommandOutput, Pos, 1) = " "
                        End If
                    Loop
                    Form1.SSPanel1.Caption = Right(GetCommandOutput, 50)
                    Call UpdateF2Prog
                End If
                If AbortFlag = 1 Then Exit Do
            End If
            
        Loop

        ' When the process terminates successfully, Err.LastDllError will be
        ' ERROR_BROKEN_PIPE (109). Other values indicates an error.

        Call CloseHandle(pi.hProcess)
    Else
        GetCommandOutput = "Failed to create process, check the path of the command line."
    End If

    ' clean up
    Call CloseHandle(hPipeRead)
    If hPipeWrite1 Then
        Call CloseHandle(hPipeWrite1)
    End If
    If hPipeWrite2 Then
        Call CloseHandle(hPipeWrite2)
    End If
    DoingShellFlag = DoingShellFlag - 1
    
End Function
