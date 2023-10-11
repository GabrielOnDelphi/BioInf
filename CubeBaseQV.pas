UNIT CubeBaseQv;

{=============================================================================================================
 Gabriel Moraru
 2016.11
 !Proprietary algorithm - Cannot be used without explicit permission!
==============================================================================================================

 Object capabilities:
    + Trim engine
    + Calculate QVs using proprietary algorithm
    + Recall N peaks

  Workflow:
  =========
   TCubeImport.Boot.ReasignPeaks;                                                                  // Fixes errors
   TCubeImport.Boot.RecallNBases;
   TCubeImport.Boot.CalculateQVs;                                                                  // Recomputes QVs for all bases

=============================================================================================================}

INTERFACE

USES
   System.SysUtils, ccCore, ccINIFile, CubeBase, CubicDNA, clRamLog, ccRichLog;

TYPE
 RTempQvData = record
    PeakSp     : Single;
    UnCaR      : Single;                                                                           { Average of UnCaR7 and UnCaR3 }
    UnCaR7     : Single;                                                                           // Uncalled/Called 7 Peak Ratio
    UnCaR3     : Single;                                                                           // Uncalled/Called 3 Peak Ratio
    PeakHeight : Single;
    Background : Single;
    Amplitude  : Single;
  end;
 TQvData= array of RTempQvData;

TYPE
 TCubeAbstractQV = class(TCubeAbstract)
  private
   FAverageDomeArea: Integer;
   FQvComputedInt: Boolean;
   function  getLowestCalledPeak    (iCurBase, NeighboursL, NeighboursR: Integer): Integer;
   function  getHighestUncalledPeak (iCurBase, NeighboursL, NeighboursR: Integer): Integer;
   function  getAverageDomeArea: Integer;                                                          { Find highest uncalled peak }
   function  recallNBase(CONST BasePos: Integer): TBase;                                           { If a new bases is detected instead of N then it returns the new base, else it returns the same N }
   function  computeDome_BaseN_NoTrace  (CONST BasePos: Integer; OUT DomeStart, DomeEnd: Integer): TBase;
  protected
   FQVExist: Boolean;
   TraceA, TraceC, TraceG, TraceT: RTrace;
   procedure FreeTraces;                                                                           { Free RAM }
   procedure ExtractTraces;
   function  traceNo2Trace(TraceNumber: Integer): RTrace;
   function  base2Trace(Base: TBase): RTrace;
   {QV}
   procedure QV_PeakSpacing  (QvData: TQvData);
   procedure QV_UnCaRatio    (QvData: TQvData);   { Uncalled peak ratio }
   procedure QV_PeakHeight   (QvData: TQvData);   { Height of this peak compared to the average peak height (on all chromatogram). (This is my own algorithm) }
   procedure QV_Background   (QvData: TQvData);   { Penalize bases that have high background }
   procedure QV_Homopolymers (QvData: TQvData);   { QHomopolymers (Multi peak separation) }                                               { How well is separated this peak from its neighbours (of the smae color). For an ideal peak its margins hould go to absolute zero. So, if two peaks are not separated by a close-to-zero region they are penalized }
   procedure MakeQVFromProbability;                                                                { The SCF does not have a single QV value. Instead it keeps 4 fields for each base (see TBase_v30). I consider that the highest filed is the final QV }
   {DOME}
   function  DomeThreshold: Real;                                                                  { If a peak dome is under x% of AverageDomeArea then ignore it (consider it background) }
   function  AverageDomeHeight  (CONST DomeStart, DomeEnd: Integer; Trace: RTrace): Cardinal;
   function  MaxDomeHeight      (CONST DomeStart, DomeEnd: Integer; Trace: RTrace): Cardinal; overload;   { Highest peak in the specified area }
   function  MaxDomeHeight      (CONST Sample: Integer;             Trace: RTrace): Cardinal; overload;
   function  ComputeAverageDomeArea: Integer;                                                      { Computes the dome area for each base then makes an average }
   function  ComputeDomeArea    (CONST Sample : Integer; Trace: RTrace; OUT DomeStart, DomeEnd: Integer): Integer; overload;
   function  ComputeDomeArea    (CONST BasePos: Integer): Integer;      overload;                  { Same as above but the input parameter is a 'base pos' not a 'sample pos' }
   function  ComputeDomeStart   (CONST Sample : Integer; Trace: RTrace): Integer;                  { Find where the 'dome' (with the center at Sample sample) starts/ends }  { Before I call this, I need to call ReasignPeaks to recompute dome center because ABI files are badly written }
   function  ComputeDomeEnd     (CONST Sample : Integer; Trace: RTrace): Integer;
   function  ComputeDome_N      (CONST Sample : Integer; OUT DomeStart, DomeEnd: Integer): TBase;  { Compute dome's start/end but for N bases. N bases are special }
   function  ComputeDome_BaseN  (CONST BasePos: Integer; OUT DomeStart, DomeEnd: Integer): TBase;
  public
   SmallDomeThreshold : Integer;                                                                                                 { In percents. If a peak dome is unde x% of AverageDomeArea then ignore it (consider it background) }
   dbgShowPeakSp      : Boolean;                                                                   { Debug. In Release version must be FALSE }
   dbgShowUnCaR       : Boolean;                                                                   { Debug. In Release version must be FALSE }
   dbgShowPeakHeight  : Boolean;                                                                   { Debug. In Release version must be FALSE }
   dbgShowBackground  : Boolean;                                                                   { Debug. In Release version must be FALSE }
   dbgShowAmplitude   : Boolean;                                                                   { Debug. In Release version must be FALSE }
   dbgShowCombined    : Boolean;                                                                   { Debug. In Release version must be TRUE  }
   EngTrim1           : REngTrim;                                                                  { the left end - for the moment it is used for both ends }
   EngTrim2           : REngTrim;                                                                  { right end }
   constructor Create(aLog: TRamLog);
   procedure Clear; override;
   {QV}
   procedure CalculateQVs;                                                                        { The procedure exits silently if QV is already calculated }
   function  AveragePeakHeight: Integer;                                                          { Media peak-urilor pt toata Chromatograma }
   function  AverageQV        : Integer;                                                          { Media QV-urilor pt bucata buna de Chromatograma }
   function  AverageQVAll     : Integer;                                                          { Arata care este media QV-urilor pt toata Chromatograma (trimmed ends included) }
   property  QVExist          : Boolean read FQVExist;                                            { True if the original chromatogram has QV data }

   function AveragePeakHeightA(Start, Stop: Integer): Integer;
   function MaxPeakHeightA: Integer;
   function NormalizeA: Integer;

   {Base caller}
   procedure RecallNBases;                                                                        {FINITE BASE CALLER}
   {OTHERS}
   procedure TrimEnds;
   procedure ReasignPeaks;                                                                        { For some strange reasons, in ABI files, the pointers between the base and the peak are placed few points BEFORE the top of the dome (the real peak). So I have to recalculate their position }
   function  Noise2SignalPercent (CONST BasePos: Integer): Single;
   property  QvComputedInternally: Boolean  read  FQvComputedInt;                                 { If true, it means that the internal algorithm was used to compute the QVs }

   property  AverageDomeArea     : Integer  read  getAverageDomeArea;
   procedure ComputeDomeMargins (CONST BasePos: Integer; OUT DomeStart, DomeEnd: Integer);        { Used by TChromaDisplay.setHilightBase }
 end;



IMPLEMENTATION
USES ccBinary, cmMath, Math;



{===============================================================================
   CREATE
===============================================================================}
constructor TCubeAbstractQV.Create;                                                               { Unless you are careful, the object might not be fully constructed when the method is called. To avoid any problems, you should override the AfterConstruction method and use that for any code that needs to wait until the object is fully constructed. If you override AfterConstruction, be sure to call the inherited method, Top. }
begin
 inherited Create(aLog);

 FQVExist:= FALSE;
 FQvComputedInt     := FALSE;
 dbgShowPeakSp      := FALSE;                                                                      { For debugging. }
 dbgShowUnCaR       := FALSE;                                                                      { For debugging. }
 dbgShowPeakHeight  := FALSE;                                                                      { For debugging. }
 dbgShowBackground  := FALSE;                                                                      { For debugging. }
 dbgShowAmplitude   := FALSE;                                                                      { For debugging. In Release version must be FALSE }
 dbgShowCombined    := TRUE;                                                                       { For debugging. In Release version must be TRUE  }
 SmallDomeThreshold := 15;                                                                         { In percents. If a peak dome is unde x% of AverageDomeArea then ignore it (consider it background) }

 InitTrimEngine(EngTrim1);                                                                         { DEFAULT TRIM ENGINE }
 InitTrimEngine(EngTrim2);
end;


procedure TCubeAbstractQV.Clear;
begin
 inherited Clear;                                                                                  { I have to call inherited. When I call clear in the code, I must also clear the parents of this object }
 FQVExist:= FALSE;
 FQvComputedInt:= FALSE;                                                                           { If true, it means that the internal algorithm was used to compute the QVs }
end;






{--------------------------------------------------------------------------------------------------
   RECOMPUTE N BASES

     Try to find out to which trace this ambiguous base belongs.
     If the trace was detected, it returns its color (base), else it returns the  original value (N).

     Applies only to N bases.
     Applies to the entire sequence (including to the gray areas).
     It is called by TCubeImport.Boot and also by TCubeAbstractSnp.DetectDoublePeaks
--------------------------------------------------------------------------------------------------}
procedure TCubeAbstractQV.RecallNBases;
VAR iCurBase: Integer;
begin
 ExtractTraces;   { Extracting the 4 traces when I need them will be too slow. Instead I precalculate them and keep them in 4 local variables }
 TRY
  for iCurBase:= 1 to NoOfBases DO
    if CharInSet(Base[iCurBase], DNA_N)
    then Base[iCurBase]:= RecallNBase(iCurBase);
 FINALLY
  FreeTraces;
 END;
end;



function TCubeAbstractQV.recallNBase (CONST BasePos: Integer): TBase;    { This is actually a base caller }
{ IMPORTANT!!! We need to call ExtractTraces before calling it! }
VAR
   ProbRec: RBaseProb;
   DomeStart, DomeEnd: Integer;
begin
 Result:= ComputeDome_BaseN_NoTrace(BasePos, DomeStart, DomeEnd);  { Guess trace associated with this N base }

 if NOT CharInSet(Result, DNA_N) then
  begin
   { Fake data for probabilities }                                  { SCF format has this feature but ABI format don't }
   FillChar(ProbRec, SizeOf(ProbRec), #0);
   case Result of
    'A': ProbRec.A:= 1;
    'C': ProbRec.C:= 1;
    'G': ProbRec.G:= 1;
    'T': ProbRec.T:= 1;
    else ProbRec.A:= 1;    {TODO 4: finish this: This is stuppid but I have to do it. I make sure that this base is assigned to a trace. }
   end;

   Assert((ProbRec.A <= 100) AND (ProbRec.c <= 100) AND (ProbRec.g <= 100) AND (ProbRec.t <= 100), 'QV out of range!');
   CellsMX[BasePos].CellProb:= ProbRec;
  end;
end;




function TCubeAbstractQV.computeDome_BaseN_NoTrace(CONST BasePos: Integer; OUT DomeStart, DomeEnd: Integer): TBase;
{------------------------------------------------------------
  IMPORTANT! I need to call ExtractTraces before calling it!
  N bases are special because in ABI files they are not assigned to a specific trace. So i have to detect which trace has the biggest dome. That dome will corespond to the true base that is hidden behind that N base. So this function returns the real base
-------------------------------------------------------------}
VAR
   Trace: RTrace;
   Average, DomeArea, CurTrace, Sample: Integer;
   DomeHeight, MaxDomeHght, MaxScore, aDomeStart, aDomeEnd: Integer;
begin
 MaxScore:= 0;                                                                                     { Biggest area for a certain 'color' means that the N base belongs to that 'color' }
 Result:= Base[BasePos];
 Sample:= Base2Smpl[BasePos];

 Assert(BasePos>= ctCellsIndex);
 Assert(BasePos<= NoOfBases);
 Assert(CharInSet(Result,DNA_N));
 Assert(traceA.Height<> NIL);

 { Find biggest dome area among all 4 traces }
 for CurTrace:= 1 to 4 DO
  begin
   Trace:= TraceNo2Trace(CurTrace);

   { Compute dome area }
   DomeArea:= ComputeDomeArea(Sample, Trace, aDomeStart, aDomeEnd);
   if DomeArea <= DomeThreshold                                                                    { If a peak dome is under x% of AverageDomeArea then ignore it (consider it background) and move to the next trace color/dome }
   then Continue;

   { Compute dome hight }
   DomeHeight := RoundEx(AverageDomeHeight(aDomeStart, aDomeEnd, Trace));
   MaxDomeHght:= MaxDomeHeight(aDomeStart, aDomeEnd, Trace);

   { Make average }
   Average:= (DomeHeight + MaxDomeHght+ DomeArea) DIV 3;                                           { The score is the average between the dome height and dome area }

   { Find the biggest dome }
   if Average> MaxScore then
     begin
      MaxScore := Average;
      Result   := Trace.Color;                                                                     { Return the real base }
      DomeStart:= aDomeStart;
      DomeEnd  := aDomeEnd;
     end;
  end;
end;


function TCubeAbstractQV.ComputeDome_N(CONST Sample: Integer; OUT DomeStart, DomeEnd: Integer): TBase; { N bases are special because in ABI files they are not assigned to a specific trace. So i have to detect which trace has the biggest dome. That dome will ocreecpond to the true base that is hidden behind that N base. So this function returns the real base }
begin
 ExtractTraces;                                                                                    { Extracting the 4 traces when I need them will be too slow. Instead I precalculate them and keep them in 4 local variables }
 TRY
   Result:= ComputeDome_BaseN_NoTrace(Sample, DomeStart, DomeEnd);
 FINALLY
  FreeTraces;
 END;
end;


function TCubeAbstractQV.ComputeDome_BaseN(CONST BasePos: Integer; OUT DomeStart, DomeEnd: Integer): TBase;
begin
 Result:= ComputeDome_N(Base2Smpl[BasePos], DomeStart, DomeEnd);
end;












{ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   CALCULATE QVs                                                                                   Documentation: http://www.insilicase.co.uk/Web/PhredScores.aspx
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - }
procedure TCubeAbstractQV.CalculateQVs;
VAR
   QvData: TQvData;
   CurBase: Integer;
   P: Single;
begin
 if NOT HasChroma
 then RAISE Exception.Create('Non chromatogram file!') at @TCubeAbstractQV.CalculateQVs;           { I need a chromatogram in order to calculate the QVs }
 if (NoOfBases< 7) then EXIT;                                                                      { Fara sta se crasuieste in QV_PeakSpacing procedure pt secvente scurte. A se vedea secventa 'amoA pCC1 220.scf' }

 { Init }
 SetLength(QvData, NoOfBases+1);                                                                   { We use a temporary matrix to store QV data }
 ExtractTraces;                                                                                    { Extracting the 4 traces when I need them will be too slow. Instead I precalculate them and keep them in 4 local variables }
 TRY
   { Compute QVs }                                                                                 { We use 5 different algorithms to compute QVs }
   QV_PeakSpacing  (QvData);                                                                       { Calculate highest and lowest peak distance in the current 7 bases window.  For each 7-bases window, the largest and smallest distance between two peaks is found and used to calculate the ratio of the largest space divided by smallest space (S(l)/S(s)).  If the sequence is evenly spaced this ratio equals 1 (the highest number is 1) }
   QV_UnCaRatio    (QvData);                                                                       { Uncalled peak ratio }
   QV_PeakHeight   (QvData);                                                                       { Height of this peak compared to the average peak height (on all chromatogram). (This is my own algorithm) }
   QV_Background   (QvData);                                                                       { Penalize bases that have high background }
   QV_Homopolymers (QvData);                                                                       { QHomopolymers (Multi peak separation) -  How well is separated this peak from its neighbours (of the smae color). For an ideal peak its margins hould go to absolute zero. So, if two peaks are not separated by a close-to-zero region they are penalized }
 FINALLY
   FreeTraces;                                                                                     { Free RAM }
 END;

 for CurBase:= ctCellsIndex to NoOfBases DO
  begin
   { All N bases get a nul QV }
   if CharInSet(Base[CurBase], DNA_N) then
    begin
     CellQV[CurBase]:= 0;
     Continue;
    end;

   if dbgShowPeakSp       then p:= QvData[CurBase].PeakSp else                                     { This is for debugging only }
   if dbgShowUnCaR        then p:= QvData[CurBase].UnCaR else                                      { This is for debugging only }
   if dbgShowPeakHeight   then p:= QvData[CurBase].PeakHeight else                                 { This is for debugging only }
   if dbgShowBackground   then p:= QvData[CurBase].Background else                                 { This is for debugging only }
   if dbgShowAmplitude    then p:= QvData[CurBase].Amplitude else                                  { This is for debugging only }
   if NOT dbgShowCombined then p:= -777                                                            { This is for debugging only }
   else
     begin
       if CurBase=13
       then EmptyDummy;

       { Combine all algorithms }                                                                  { The final QV score is the worst number returned by any of the 5 algorithms }
       p:= Min3S(QvData[CurBase].Amplitude, QvData[CurBase].Background, QvData[CurBase].PeakHeight);
       p:= Min3S(p, QvData[CurBase].UnCaR, QvData[CurBase].PeakSp);
     end;

   { Finally we have the QV }
   p:= Ensure100(p);
   CellQV[CurBase]:= RoundEx(p);
  end;

 SetLength(QvData, 0);                                                                             { Release memory }
 FQVExist:= TRUE;
 FQvComputedInt:= TRUE;                                                                            { If true, it means that the internal algorithm was used to compute the QVs }
 RamLog.AddVerb('  QVs computed with DNA Baser''s internal algorithm.');
end;



{ QV STEP 1: }
{ QHomopolymers (Multi peak separation) }                                                          { How well is separated this peak from its neighbours (of the smae color). For an ideal peak its margins hould go to absolute zero. So, if two peaks are not separated by a close-to-zero region they are penalized }
procedure TCubeAbstractQV.QV_Homopolymers;
VAR
   CurTrace: RTrace;
   CurBase: TBase;
   Center, Dif, MinAmpl, RightAmpl, LeftAmpl, LeftDome, RightDome, iCurBase, MainPeakHeight: Integer;
begin
 { Peak identification for the four traces }
 for iCurBase:= ctCellsIndex to NoOfBases DO
  begin
   CurBase:= UpCase(Base[iCurBase]);

   { Ambiguity? }
   if NOT CharInSet(CurBase, DNA_ACGT) then
    begin
     QvData[iCurBase].Amplitude:= 1;
     Continue;
    end;

   { Obtain current trace }
   CurTrace:= Base2Trace(CurBase);

   Center:= Base2Smpl[iCurBase];
   MainPeakHeight:= CurTrace.Height[Center];
   if MainPeakHeight= 0
   then MainPeakHeight:= 1;                                                                        { To prevent 'Division by zero' }

   Assert(Center<= NoOfSamples);
   Assert(Center>= 1);

   LeftDome := ComputeDomeStart(Center, CurTrace);
   RightDome:= ComputeDomeEnd  (Center, CurTrace);

   Assert(LeftDome<= NoOfSamples);
   Assert(LeftDome>= 1);

   LeftAmpl := CurTrace.Height[LeftDome];
   RightAmpl:= CurTrace.Height[RightDome];

   { Find the worst neighbor }                                                                     { This is the part of the dome between the main trace and its neighbors that refuses to go down to zero }
   MinAmpl:= max(LeftAmpl, RightAmpl);
   Dif:= MainPeakHeight - MinAmpl;                                                                 { The difference between the height of the main peak and its worst neighbor }

   if Dif< 0
   then QvData[iCurBase].Amplitude:= 0
   else QvData[iCurBase].Amplitude:= (Dif * 100) / MainPeakHeight;
  end;
end;



{ QV STEP 2: Signal/noise ratio }
procedure TCubeAbstractQV.QV_Background;                                                           { Penalize bases that have high background }
VAR
   CurBase: TBase;
   iCurBase: Integer;
   BkgNoiseProc: Single;
begin
 for iCurBase:= ctCellsIndex to NoOfBases DO
  begin
   CurBase:= UpCase(Base[iCurBase]);

   if NOT CharInSet(CurBase,DNA_ACGT) then
    begin
     QvData[iCurBase].Background:= 0;      {TODO 5: finish this: I have found a polimorphic (IUPAC) base }  { Peak identification for the four traces }
     Continue;
    end;

   BkgNoiseProc:= Noise2SignalPercent(iCurBase);

   { Scale it down }
   if BkgNoiseProc<= 0 then BkgNoiseProc:= 0.1;                                                   { The LN function does not accept zero as input }
   BkgNoiseProc:= 0.0059*BkgNoiseProc*BkgNoiseProc -1.5537*BkgNoiseProc+ 71;                      { original formula: 0.0059*BkgNoiseProc*BkgNoiseProc -1.5537*BkgNoiseProc+ 97.287 }
   if BkgNoiseProc> 100 then BkgNoiseProc:= 100;

   QvData[iCurBase].Background:= BkgNoiseProc;
  end;
end;


function TCubeAbstractQV.Noise2SignalPercent(CONST BasePos: Integer): Single;                     { The percent of bck area that overlaps (under) the signal }
VAR
   Sample, Smpl, BackgroundArea: Integer;
   CurBase: TBase;
   RefBaseArea, MaxTrace, iStart, iEnd: Integer;
begin
 Result:= 0;
 Assert(BasePos>= ctCellsIndex);
 Assert(BasePos<= noofbases);

 { Peak identification for the four traces }
 CurBase:= UpCase(Base[BasePos]);

 { Prepare 'N' bases }                                                                            { I have found a polimorphic (IUPAC) base. THIS SHOULD NEVER HAPPEN. The caller should take care of it. }
 if NOT CharInSet(CurBase, DNA_ACGT) then EXIT;

 Sample:= Base2Smpl[BasePos];

 RefBaseArea:= 0;
 BackgroundArea:= 0;
 case CurBase of
  'A': begin
        RefBaseArea:= ComputeDomeArea(Sample, TraceA, iStart, iEnd);                              { Compute area of main peak }
        for Smpl:= iStart to iEnd DO
         begin
          MaxTrace:= Find_Max(TraceC.Height[Smpl], TraceG.Height[Smpl], TraceT.Height[Smpl]);
          if MaxTrace> TraceA.Height[smpl]
          then BackgroundArea:= BackgroundArea+ TraceA.Height[smpl]
          else BackgroundArea:= BackgroundArea+ MaxTrace;
        end;
       end;

  'C': begin
        RefBaseArea:= ComputeDomeArea(Sample, TraceC, iStart, iEnd);
        for Smpl:= iStart to iEnd DO
         begin
          MaxTrace:= Find_Max(TraceT.Height[Smpl], TraceG.Height[Smpl], TraceA.Height[Smpl]);
          if MaxTrace> TraceC.Height[smpl]
          then BackgroundArea:= BackgroundArea+ TraceC.Height[smpl]
          else BackgroundArea:= BackgroundArea+ MaxTrace;
        end;
       end;

  'G': begin
        RefBaseArea:= ComputeDomeArea(Sample, TraceG, iStart, iEnd);
        for Smpl:= iStart to iEnd DO
         begin
          MaxTrace:= Find_Max(TraceC.Height[Smpl], TraceA.Height[Smpl], TraceT.Height[Smpl]);
          if MaxTrace> TraceG.Height[smpl]
          then BackgroundArea:= BackgroundArea+ TraceG.Height[smpl]
          else BackgroundArea:= BackgroundArea+ MaxTrace;
         end;
        end;

  'T': begin
        RefBaseArea:= ComputeDomeArea(Sample, TraceT, iStart, iEnd);
        for Smpl:= iStart to iEnd DO
         begin
          MaxTrace:= Find_Max(TraceC.Height[Smpl], TraceG.Height[Smpl], TraceA.Height[Smpl]);

          if MaxTrace> TraceT.Height[smpl]
          then BackgroundArea:= BackgroundArea+ TraceT.Height[smpl]
          else BackgroundArea:= BackgroundArea+ MaxTrace;
        end;
       end;
 end;

 if RefBaseArea= 0
 then Result:= 0
 else Result:= (BackgroundArea / RefBaseArea)* 100;
end;



{ QV STEP 3: Peak height }                                          { Height of this peak compared to the average peak height (on all chromatogram). (This is my own algorithm) }
procedure TCubeAbstractQV.QV_PeakHeight;
VAR
   AvPeakHeight, CurBase: Integer;
   PeakHeight: Double;
begin
 AvPeakHeight:= AveragePeakHeight;

 if AvPeakHeight <= 0                                               { This happens for really bad/short samples }
 then
  for CurBase:= ctCellsIndex to NoOfBases
   DO QvData[CurBase].PeakHeight:= 0

 else
  for CurBase:= ctCellsIndex to NoOfBases do
   begin
    PeakHeight:= Base2SampleHeight(CurBase);
    PeakHeight:= ProcentRepresent(PeakHeight, AvPeakHeight);        { div 2 to scale down }

    if PeakHeight> 100
    then PeakHeight:= 100;

    QvData[CurBase].PeakHeight:= PeakHeight;
   end;
end;




{ QV STEP 4: Uncalled to called peak height ratio }                 { For each x-base window the ratio of the height of the largest uncalled peak (P(lu)) against the height of the shortest called peak (P(sc)) is found. If no uncalled peak is present in the window then the highest background value under a peak is used    Explanation source 2:  Window: 3/7 peaks around the current one. Calculation: (amplitude of the largest uncalled peak) / (smallest called peak).An uncalled peak is a peak in the signal that was not assigned to a predicted location by phred and thus does not result in a base call. If the called base is an N, phred assigns a poor value.  Note that this is not what is sometimes called the signal to noise ratio, as uncalled peaks may be true peaks missed by the base-calling program rather than noise in the conventional sense. The minimum parameter value is 0 for traces with no uncalled peaks.
{ SUB1. Find lowest called peak }
function TCubeAbstractQV.getLowestCalledPeak(iCurBase, NeighboursL, NeighboursR: Integer): Integer;
VAR
   CurBaseInWnd, CurSmpl: Integer;
   aBase: TBase;
begin
 Result:= 1;                                                         { It will return 1 for non ACGT bases }

 for CurBaseInWnd:= (iCurBase- NeighboursL) to (iCurBase+ NeighboursR) DO
  begin
   aBase:= UpCase(Base[CurBaseInWnd]);
   CurSmpl:= Base2Smpl[CurBaseInWnd];

   if (aBase= 'A') AND (Result < Sample[CurSmpl].HeightA)
   then Result:= Sample[CurSmpl].HeightA else

   if (aBase= 'G') AND (Result < Sample[CurSmpl].HeightG)
   then Result:= Sample[CurSmpl].HeightG else

   if (aBase= 'C') AND (Result < Sample[CurSmpl].HeightC)
   then Result:= Sample[CurSmpl].HeightC else

   if (aBase= 'T') AND (Result < Sample[CurSmpl].HeightT)
   then Result:= Sample[CurSmpl].HeightT;
  end;

 if Result= 0
 then Result:= 1;
end;

{ SUB2. Find highest uncalled peak }                                          { Compare the height of the trace associated with iCurBase base with all other 3 traces }
function TCubeAbstractQV.getHighestUncalledPeak(iCurBase, NeighboursL, NeighboursR: Integer): Integer;
VAR
  CurBase: TBase;
  SmplHeight, MainSmplHeight, CurSmpl, SmplStart, SmplEnd: Integer;
begin
 Result:= 0;
 CurBase:= UpCase(Base[iCurBase]);

 { Calculate the 3/7 base window }
 SmplStart:= Base2Smpl[iCurBase- NeighboursL];
 SmplEnd  := Base2Smpl[iCurBase+ NeighboursR];

 { A }
 if CurBase<> 'A' then                                                        { Do not compare A with A, C with C and so on... }
  for CurSmpl:= SmplStart+2 to SmplEnd-2 DO                                   { Compare the height of this base with trace A }
   begin
      SmplHeight:= Sample[CurSmpl].HeightA;
      MainSmplHeight:= Base2SampleHeight(iCurBase, CurSmpl);                  { Sample height belonging to the base which we compute now (iBasePos) }
      if (SmplHeight <= MainSmplHeight) AND (SmplHeight >  Result)            { If this trace (Y) is under the trace associated to the 'iCurBase' base, we store it as the 'highest uncalled peak' }
      then Result:= SmplHeight;
   end;

 { C }
 if CurBase<> 'C' then                                                        { Do not compare A with A, C with C and so on... }
  for CurSmpl:= SmplStart+2 to SmplEnd-2 DO                                   { Compare the height of this base with trace C }
   begin
      SmplHeight:= Sample[CurSmpl].HeightC;
      MainSmplHeight:= Base2SampleHeight(iCurBase, CurSmpl);
      if (SmplHeight > Result) AND (SmplHeight <= MainSmplHeight)
      then Result:= SmplHeight;
   end;

 { G }
 if CurBase<> 'G' then
  for CurSmpl:= SmplStart+2 to SmplEnd-2 DO
   begin
      SmplHeight:= Sample[CurSmpl].HeightG;
      MainSmplHeight:= Base2SampleHeight(iCurBase, CurSmpl);
      if (SmplHeight > Result) AND (SmplHeight <= MainSmplHeight)
      then Result:= SmplHeight;
   end;

 { T }
 if CurBase<> 'T' then
  for CurSmpl:= SmplStart+2 to SmplEnd-2 DO
   begin
      SmplHeight:= Sample[CurSmpl].HeightT;
      MainSmplHeight:= Base2SampleHeight(iCurBase, CurSmpl);
      if (SmplHeight > Result) AND (SmplHeight<= MainSmplHeight)
      then Result:= SmplHeight;
   end;

 if Result= 0
 then Result:= 1;                                                                                  { Prevent 'Division by zero' }
end;


procedure TCubeAbstractQV.QV_UnCaRatio;  { Uncalled peak ratio }
CONST
   ctMaxNeighbors7 = 3;
VAR
   Value: Single;
   CurBase, NeighboursL, NeighboursR, LowestCalledPeak, HighUncalledPeak: Integer;
begin
 {TODO 5: If there is no uncalled peak, the largest of the three uncalled trace values at the location of the called base peak is used instead. }

 { 7-base window }
 for CurBase:= ctCellsIndex to NoOfBases DO
  begin
   { Precalculate right and left neighbors }
   if (CurBase > 3) AND (CurBase < NoOfBases- 3)                                                   { Check if we have 3 neighbors on each side}
   then
    begin
     NeighboursL := ctMaxNeighbors7;
     NeighboursR := ctMaxNeighbors7;
    end
   else
     if CurBase <= ctMaxNeighbors7                                                                 { Check if the base is at the beginning (left) of the array }
     then
       begin
        NeighboursL := CurBase-1;
        NeighboursR := 3;
       end
     else                                                                                          { The base is at the right end of the array }
       begin
        NeighboursR := NoOfBases- CurBase;
        NeighboursL := 3;
       end;

   if CharInSet(Base[CurBase],DNA_N)
   then QvData[CurBase].UnCaR7:= 0
   else
    begin
     { Calculate ratio }
     LowestCalledPeak := getLowestCalledPeak(CurBase, NeighboursL, NeighboursR);                   { Find lowest called peak }
     HighUncalledPeak := getHighestUncalledPeak(CurBase , NeighboursL, NeighboursR);               { Find highest uncalled peak }
     QvData[CurBase].UnCaR7:= 100 - (HighUncalledPeak/LowestCalledPeak) * 100;                     { Small values are ideal values }
    end;
  end;

 { 3-base window }
 for CurBase:= ctCellsIndex to NoOfBases DO
  begin
   { Precalculate right and left neighbors }
   if (CurBase > 1) AND (CurBase < NoOfBases- 1)                                                   { Check if we have 1 neighbors on each side}
   then
    begin
     NeighboursL := 1;
     NeighboursR := 1;
    end
   else
     if CurBase <= 1                                                                               { Check if the base is at the beginning (left) of the array }
     then
       begin
        NeighboursL := CurBase-1;
        NeighboursR := 1;
       end
     else                                                                                          { The base is at the right end of the array }
       begin
        NeighboursR := NoOfBases- CurBase;
        NeighboursL := 1;
       end;

   { Calculate ratio }
   if CharInSet(Base[CurBase], DNA_N)
   then QvData[CurBase].UnCaR3:= 0
   else
    begin
     LowestCalledPeak := GetLowestCalledPeak   (CurBase, NeighboursL, NeighboursR);                { Find lowest called peak }
     HighUncalledPeak := GetHighestUncalledPeak(CurBase, NeighboursL, NeighboursR);                { Find highest uncalled peak }
     QvData[CurBase].UnCaR3:= 100 - (HighUncalledPeak/LowestCalledPeak) * 100;
    end;
  end;

 { Make average between QvData[CurBase].UnCaR3 and QvData[CurBase].UnCaR7 }
 for CurBase:= 1 to NoOfBases DO
  begin
   Value:= (QvData[CurBase].UnCaR7 + QvData[CurBase].UnCaR3) / 2;
   if Value> 100
   then Value:= 100;
   QvData[CurBase].UnCaR:= Value;
  end;
end;




{ QV STEP 5. Peak spacing }    { Calculate highest and lowest peak distance in the current 7 bases window.  For each 7-bases window, the largest and smallest distance between two peaks is found and used to calculate the ratio of the largest space divided by smallest space (S(l)/S(s)).  If the sequence is evenly spaced this ratio equals 1 (the highest number is 1) }
procedure TCubeAbstractQV.QV_PeakSpacing(QvData: TQvData);
CONST ctMaxNeighbors = 3;
VAR
   NeighboursL, NeighboursR, smpl1, smpl2, dist: Integer;
   CurBaseInWnd, CurBase, Max_Peak_Space, Min_Peak_Space: Integer;
begin
 for CurBase:= ctCellsIndex to NoOfBases DO
  begin
   { Precalculate right and left neighbors }
   if (CurBase > ctMaxNeighbors) AND (CurBase < NoOfBases- ctMaxNeighbors)                         { Check if we have 3 neighbors on each side}
   then
    begin
     NeighboursL := ctMaxNeighbors;
     NeighboursR := ctMaxNeighbors;
    end
   else
     if CurBase <= ctMaxNeighbors                                                                  { Check if the base is at the beginning of the array}
     then
       begin
        NeighboursL := CurBase-1;
        NeighboursR := ctMaxNeighbors;
       end
     else                                                                                          { The base is at the right end of the array }
       begin
        NeighboursR := NoOfBases- CurBase;
        NeighboursL := ctMaxNeighbors;
       end;

   Max_Peak_Space := 0;
   Min_Peak_Space := High(Integer);

   for CurBaseInWnd:= (CurBase- NeighboursL) to (CurBase+ NeighboursR)-1 DO                        { Apply the algorithm in the current 7-base window }
    begin
     { Find distance between these 2 peaks }
     smpl1 := Base2Smpl[CurBaseInWnd  ];
     smpl2 := Base2Smpl[CurBaseInWnd+1];
     dist  := smpl2- smpl1;

     if Max_Peak_Space< dist
     then Max_Peak_Space:= dist;
     if Min_Peak_Space> dist
     then Min_Peak_Space:= dist;
    end;

   if Max_Peak_Space<= 0
   then Max_Peak_Space:= 1;                                                                        { For some samples I get 'invalid floating point operation' so I have to make sure I don't get it }

   { Finally assign the new calculated QV }
   QvData[CurBase].PeakSp:= (Min_Peak_Space/Max_Peak_Space) * 100;                                 { The numbers are in 0-1 range. So I multiply with 100 to normalize them. }
 end;
end;





{--------------------------------------------------------------------------------------------------
   REASIGN PEAKS
--------------------------------------------------------------------------------------------------}
procedure TCubeAbstractQV.ReasignPeaks;                                                            { For some strange reasons, in ABI files, the pointers between the base and the peak are placed few points BEFORE the top of the dome (the real peak). So I have to recalculate their position }
CONST
   MaxDistance= 5;
VAR
   aBase: TBase;
   Ascending, Descending: Boolean;
   Left, Center, Right, iTo, iDownTo, iHeight, iHeight2, OldPointer, CurBase, CurSmpl: Integer;
begin
 for CurBase:= ctCellsIndex to NoOfBases DO
  begin
   aBase:= Base[CurBase];
   if NOT CharInSet(abase, DNA_ACGT) then Continue;                                                { The N/IUPAC bases are not assigned to a certain trace so they cannot be reasigned }

   OldPointer:= Base2Smpl[CurBase];

   { Make sure I don't go out of range }
   if OldPointer-1< ctChromaIndex then Continue;
   if OldPointer+1> NoOfSamples   then Continue;

   { On which slope am I? }

   Left  := Base2SampleHeight (aBase, OldPointer-1);
   Center:= Base2SampleHeight (aBase, OldPointer  );
   Right := Base2SampleHeight (aBase, OldPointer+1);


   Ascending := (Left< Center) AND (Center< Right);
   Descending:= (Left> Center) AND (Center> Right);

   { I am already,a peak or I am on a plateau }
   if Ascending = Descending then Continue;

   if Ascending                                                                                    { Search for the peak from this point to the right }
   then
    begin
     iTo:= OldPointer+ MaxDistance;
     if iTo> NoOfSamples -1 then EXIT;                                                             { Make sure I don't go out of range }
     for CurSmpl:= OldPointer to iTo DO
      begin
       iHeight := Base2SampleHeight(aBase, CurSmpl   );
       iHeight2:= Base2SampleHeight(aBase, CurSmpl+ 1);

       if iHeight>= iHeight2 then                                                                  { aici trebuie sa fie neaparat 'mai mare sau egal' }
        begin
          if Chroma [CurSmpl].Ptr2Base= NoneAssigned then                                          { If there is already a base assigned to this peak do nothing. I already discussed this with Cristina and agreed on that }
           begin
            Chroma [OldPointer].Ptr2Base:= NoneAssigned;                                           { Delete old pointer }
            Chroma [CurSmpl].Ptr2Base:= CurBase;                                                   { Reasign pointers }
            CellsMX[CurBase].Ptr2Smpl:= CurSmpl;
           end;
          Break;
        end;
      end;
    end
   else  { Descending }                                                                            { Search for the peak from this point to the left }
     begin
      iDownTo:= OldPointer- MaxDistance;
      if iDownTo< ctCellsIndex+ 1 then EXIT;
      for CurSmpl:= OldPointer downto iDownTo DO
       begin
        iHeight := Base2SampleHeight(aBase, CurSmpl   );
        iHeight2:= Base2SampleHeight(aBase, CurSmpl- 1);

        if iHeight>= iHeight2 then
         begin
           if Chroma [CurSmpl].Ptr2Base= NoneAssigned then
            begin
             Chroma [OldPointer].Ptr2Base:= NoneAssigned;                                          { Delete old pointer }
             Chroma [CurSmpl].Ptr2Base:= CurBase;                                                  { Reasign pointers }
             CellsMX[CurBase].Ptr2Smpl:= CurSmpl;
            end;
           Break;
         end;
       end;
     end;
  end;
end;




{--------------------------------------------------------------------------------------------------
   DOME START/END
--------------------------------------------------------------------------------
 Find where the 'dome' (with the center at Sample sample) starts/ends.
 Before I call this, I need to call ReasignPeaks to recompute dome center because ABI files are badly written.
 If Result= -1 it means that that the end of the dome cannot be found          }

function TCubeAbstractQV.ComputeDomeStart(CONST Sample: Integer; Trace: RTrace): Integer;
VAR LeftNeighbor, smpl: Integer;
begin
 Assert(Length(Trace.Height) > 0);
 Assert(Sample> 0, i2s(Sample));
 Assert(Sample<= NoOfSamples, i2s(Sample));
 Assert(CharInSet(Trace.Color,DNA_ACGTN));
 Result:= Sample;                                                                                  { This is necessary for some special cases such as DomeStart= 1 }

 { Find left neighbor I }
 for smpl:= Sample downto ctChromaIndex+1 DO                                                       { Initially, both the start and end point of the peak region are equal to the reference point (p) }
   if Trace.Height[smpl-1]>= Trace.Height[smpl] then                                               { The point is moved apart from the reference point (p), until the height at the left of the start point is higher than the height at the start point }
    begin
     Result:= smpl;
     Break;
    end;

 if Sample= ctChromaIndex then EXIT;                                                               { This is the first sample. No more neighbors on its left side }

 { Find left neighbor II }
 LeftNeighbor:= 0;
 for smpl:= Sample-1 downto ctChromaIndex DO
  if SampleHasBaseAssignedEx(smpl, Trace.Color) then
    begin
     LeftNeighbor:= smpl;
     Break;
    end;

 { Don't let it step over the left neighbor }                                                      { This can happen at the ends when a single long dome has several bases assigned to it. The bases are not separate through 'valeys' }
 if Result<= LeftNeighbor
 then Result:= LeftNeighbor+ ((LeftNeighbor- Result) DIV 2);

 if Result<= LeftNeighbor                                                                          { One extrac check to see if the Result was indeed fixed by the above operation }
 then Result:= LeftNeighbor+ 1;

 if Result> Sample
 then Result:= Sample;
end;


function TCubeAbstractQV.ComputeDomeEnd(CONST Sample: Integer; Trace: RTrace): Integer;            { If Result= -1 it means that that the end of the dome cannot be found }
VAR RightNeighbor, smpl: Integer;
begin
 Assert(Length(Trace.Height)= Length(Chroma));
 Assert(Sample>= ctChromaIndex);
 Assert(Sample<= Length(Trace.Height));

 Result:= Sample;
 for smpl:= Sample to Length(Trace.Height)- IndexedIn1- 2 DO                                       { -2 pt ca am un pic mai jos +2}
   if  (Trace.Height[smpl  ]>= Trace.Height[smpl+1])
   AND (Trace.Height[smpl+1]<= Trace.Height[smpl+2]) then
    begin
     Result:= smpl+ 1;
     Break;
    end;

 if Sample= NoOfBases then EXIT;                                                                   { No more neighbors on its left side }

 { Find right neighbor }
 RightNeighbor:= 0;
 for smpl:= Sample+1 to NoOfSamples DO
  if SampleHasBaseAssignedEx(smpl, Trace.Color) then
   begin
    RightNeighbor:= smpl;
    Break;
   end;

 { Don't let it step over the right neighbor }                                                     { This can happen at the ends when a single long dome has several bases assigned to it. The bases are not separate through 'valeys' }
 if Result>= RightNeighbor
 then Result:= RightNeighbor- ((Result- RightNeighbor) DIV 2);

 if Result>= RightNeighbor
 then Result:= RightNeighbor- 1;

 if Result< Sample
 then Result:= Sample;
end;





{--------------------------------------------------------------------------------------------------
   DOME AREA
--------------------------------------------------------------------------------------------------}
function TCubeAbstractQV.ComputeAverageDomeArea: Integer;                                          { Computes the dome area for each base then makes an average }
VAR CurBase: Integer;
begin
 Result:= 0;
 for CurBase:= ctCellsIndex to NoOfBases
  DO Result:= Result+ ComputeDomeArea(CurBase);
 Result:= Result DIV NoOfBases;
end;


function TCubeAbstractQV.ComputeDomeArea(CONST Sample: Integer; Trace: RTrace; OUT DomeStart, DomeEnd: Integer): Integer;
VAR Smpl: Integer;
begin
 { Find margins }
 DomeStart := ComputeDomeStart (Sample, Trace);
 if (DomeStart<= 0) then EXIT(0);

 DomeEnd:= ComputeDomeEnd(Sample, Trace);
 if (DomeEnd<= 0) then EXIT(0);

 { Compute area }
 Result:= 0;
 for Smpl:= DomeStart to DomeEnd DO
  Result:= Result+ Trace.Height[Smpl];
end;


function TCubeAbstractQV.ComputeDomeArea(CONST BasePos: Integer): Integer;                         { Same as above but the input parameter is a 'base pos' not a 'sample pos' }
VAR
   Sample: Integer;
   DomeStart, DomeEnd: Integer;
   Trace: RTrace;
begin
 Sample:= Base2Smpl[BasePos];
 Trace := Base2Trace(Base[BasePos]);

 if Trace.Height= NIL                                                                              { Trace is NIL for non ACGT bases because it is difficult to tell exactly to which trace that base belongs. Theoretically N bases belong to all 4 traces, and UIPAC bases belong to one or more traces. But this will complicate the code toomuch }
 then Result:= 0
 else Result:= ComputeDomeArea(Sample, Trace, DomeStart, DomeEnd)
end;


function TCubeAbstractQV.getAverageDomeArea: Integer;                                              { Property getter }
begin
 if FAverageDomeArea= 0
 then FAverageDomeArea:= ComputeAverageDomeArea;                                                   { If the AverageDomeArea was never calculated, calculate it now }
 Result:= FAverageDomeArea;
end;



{--------------------------------------------------------------------------------------------------
   DOME MARGINS
   Similar to ComputeDomeArea but it ignores the Area and returns only the margins
--------------------------------------------------------------------------------------------------}

procedure TCubeAbstractQV.ComputeDomeMargins(CONST BasePos: Integer; OUT DomeStart, DomeEnd: Integer);
VAR Trace: RTrace;
    Sample: Integer;
begin
 Assert(CharInSet(Base[BasePos],DNA_ACGT));                                                        { This is because ExtractTrace can only work on one of the 4 traces. 'N' bases are not connected to any of the 4 traces }
 Assert(BasePos> 0, i2s(BasePos));
 Assert(BasePos<= NoOfBases, i2s(BasePos));

 if CharInSet(Base[BasePos], DNA_ACGT)
 then
  begin
   ExtractTrace(Base[BasePos], Trace);                                                             { Extracting the 4 traces when I need them is slow! }

   Sample:= Base2Smpl[BasePos];
   DomeStart:= ComputeDomeStart (Sample, Trace);
   DomeEnd  := ComputeDomeEnd   (Sample, Trace);

   SetLength(Trace.Height, 0);
  end
 else ComputeDome_BaseN(BasePos, DomeStart, DomeEnd);
end;







{===============================================================================
   TRIM ENGINE
===============================================================================}
procedure TCubeAbstractQV.TrimEnds;
VAR
  i1, i2, i3, x, j, iGoodQVStart, iGoodQVEnd: Integer;
  Prag: real;
  NrOfGoodBases: Integer;                                                                          { Numarul de baze care au QV bun (peste prag) in fereastra data }
  TempQVStart  : Integer;                                                                          { Poz unde am gasit prima baza care depaseste pragul. Pentru secvente incredibil de rele, s-ar putea ca niciodata sa nu am o baza sau fereastra care atinge pragul }
begin
 Assert(EngTrim1.GoodBasesNr> 0);                                                                  { Make sure I didn't forgot to set this }
 Assert(EngTrim1.WindLength > 0);
 Assert(EngTrim1.GoodQVTresh> 0);

 { INIT }
 DirtyGoodBases:= TRUE;
 TempQVStart :=  0;
 iGoodQVStart:= -1;
 iGoodQVEnd  := -1;
 Prag:= (EngTrim1.WindLength* EngTrim1.GoodBasesNr) / 100;

 if FQVExist
 AND (EngTrim1.GoodQVTresh> 1)                                                                     { If GoodQVTresh= 1 the I consider that I don't want to apply the trimming }
 then
  BEGIN {#1}

  {----------------* CURAT BAZELE DE LA CAP *-----------------}

  { A. CAUT PRIMA BAZA BUNA }
  for i1:= ctCellsIndex to NoOfBases DO
   if CellQV[i1] >= EngTrim1.GoodQVTresh then
    begin
     TempQVStart:= i1;                                                                             { memoreaza de unde incepe prima baza buna }
     Break;
    end;

  { B. CAUT PANA CAND FEREASTRA CONTINE DESTULE BAZE BUNE }
  if TempQVStart> 0 then                                                                           { TEST AGANST EXTEME LOW QV}
   begin
    for i2:= TempQVStart TO (NoOfBases- EngTrim1.WindLength) DO
     begin
       NrOfGoodBases:= 0;

       { OBTIN CALITATEA FERESTREI }
       for j:= i2 TO (i2+EngTrim1.WindLength) DO                                                   { aflu cate baze sunt de calitate in fereastra }     { GoodQVWindowSize= CATE baze cu QV bun sa am ca sa declar fereastra OK }
          if CellQV[j] >= EngTrim1.GoodQVTresh
          then inc(NrOfGoodBases);

       if NrOfGoodBases>= Prag then                                                                { daca mai mult de x % din baze sunt bune    ->   x=80% }
        begin
         iGoodQVStart:= i2;                                                                        { memoreaza de unde incepe fereastra buna (a nu se confunda cu prima baza buna) }
         Break;
        end;
     end;
   end
  else iGoodQVStart:= -1;                                                                          { ALARMA! Nu am gasit nici o baza buna }

  { C. CAUT INAUNTRU FERESTREI PRIMA BAZA CU ADEVARAT BUNA }
  IF iGoodQVStart> 0 THEN
   BEGIN{@}                                                                                        { TEST AGANST EXTEME LOW QV}
    for i3:= iGoodQVStart TO NoOfBases DO
     if CellQV[i3] >= EngTrim1.GoodQVTresh then
      begin
       iGoodQVStart:= i3;                                                                          { memoreaza de unde incepe prima baza buna }
       break;
      end;

    {----------------* CURAT BAZELE DE LA COADA *----------------}

    { A. CAUT PRIMA BAZA BUNA }
    for x:= NoOfBases downto 1 DO
      if CellQV[x] >= EngTrim1.GoodQVTresh then
       begin
        iGoodQVEnd:= x;                                                                            { memoreaza de unde incepe prima baza buna }
        break;
       end;
    { B. CAUT PANA CAND FEREASTRA CONTINE DESTULE BAZE BUNE }
    for x:= iGoodQVEnd downto (1+EngTrim1.WindLength) do
     begin
      NrOfGoodBases:= 0;
      { OBTIN CALITATEA FERESTREI }
      for j:= x downto (x-EngTrim1.WindLength)
       DO if CellQV[j] >= EngTrim1.GoodQVTresh
          then inc(NrOfGoodBases);

      if NrOfGoodBases>= Prag  then                                                                { daca mai mult de x% baze sunt bune }
       begin
        iGoodQVEnd:= x;                                                                            { memoreaza de unde incepe prima baza buna }
        Break;
       end;
     end;
     { C. CAUT INAUNTRU FERESTREI PRIMA BAZA CU ADEVARAT BUNA }
    for x:= iGoodQVEnd downto 1 DO
     if CellQV[x] >= EngTrim1.GoodQVTresh then
      begin
       iGoodQVEnd:= x;                                                                             { memoreaza de unde incepe prima baza buna }
       break;
      end;
   END{@}
  ELSE
   begin                                                                                           { ALARMA! Nu am gasit nici o baza buna }
     iGoodQVStart:= NoOfBases-1;                                                                   { pentru cazul in care TOATE bazele sunt rele }
     iGoodQVEnd  := NoOfBases;                                                                     { pentru cazul in care TOATE bazele sunt rele }
   end;

  { CHANGED }
  { REBUILD stuff pentru ca obiectul a fost modificat }
  GoodQVStart:= iGoodQVStart;                                                                      { This will call also buildGoodBases }
  LastGoodBase:= iGoodQVEnd;                                                                       { This will call also buildGoodBases }
 END {#1}

 ELSE {QVs does not exist}
  begin
   GoodQVStart:= 1;
   LastGoodBase:= NoOfBases;
  end;
end;





{--------------------------------------------------------------------------------------------------
   QV - Average Peak Height
--------------------------------------------------------------------------------------------------}
function TCubeAbstractQV.AveragePeakHeight: Integer;                                               { Media peak-urilor (care au baze asignate) pt toata cromatograma (all 4 traces) }
VAR Sum, CurSample, CurBase, Biggest: Integer;
begin
 Assert(HasChroma);

 Sum:= 0;
 for CurBase:= ctCellsIndex to NoOfBases DO
  if Base[CurBase]<> Gap then                                                                      { Atunci cand calculez calitatea, nu numar si gap-urile introduse de asamblare }
   begin
    CurSample:= Base2Smpl[CurBase];

    { Invalid sample protection }
    if (CurSample< ctChromaIndex)                                                                  { NOTA: 7 useri au raportat assertion failure aici }
    { KEEP THE BREAKPOINT HERE }then CurSample:= ctChromaIndex;

    if (CurSample> NoOfSamples)
    { KEEP THE BREAKPOINT HERE }then CurSample:= NoOfSamples;

    Biggest:= 0;                                                                                   { aflu care din cele 4 traceuri/peakuri asociate bazei curente, are valoarea cea mai mare (e cel mai inalt) }
    if Biggest< Chroma[CurSample].HeightA then biggest:= Chroma[CurSample].HeightA;
    if Biggest< Chroma[CurSample].HeightC then biggest:= Chroma[CurSample].HeightC;
    if Biggest< Chroma[CurSample].HeightG then biggest:= Chroma[CurSample].HeightG;
    if Biggest< Chroma[CurSample].HeightT then biggest:= Chroma[CurSample].HeightT;
    Sum:= Sum+ Biggest;
   end;

 Result:= Sum div NoOfBases;
end;



function TCubeAbstractQV.MaxPeakHeightA: Integer;                                                  { Biggest peak for A trace }
VAR Smpl: Integer;
begin
 Assert(HasChroma);

 Result:= 0;
 for Smpl:= 1 to NoOfSamples DO
   begin
    if Result < Chroma[Smpl].HeightA
    then Result:= Chroma[Smpl].HeightA;
   end;
end;


function TCubeAbstractQV.AveragePeakHeightA(Start, Stop: Integer): Integer;                        { Media peak-urilor pt trace-ul A }
VAR Sum: Cardinal;
    Smpl: Integer;
begin
 Sum:= 0;

 for Smpl:= Start to Stop
  DO Sum:= Sum+ Chroma[Smpl].HeightA;

 Result:= Sum div Cardinal(Stop-Start+1);
end;


function TCubeAbstractQV.NormalizeA: Integer;                                                      { Media peak-urilor pt trace-ul A }
VAR
  Smpl: Integer;
  WndStart: Integer;
  WndEnd: Integer;
  WndAverage: Integer;
  AverageHeight: Integer;
  Mean, MaxHeight: Integer;
begin
 Result:= -77;

 AverageHeight:= AveragePeakHeight;
 MaxHeight    := MaxPeakHeightA;
 Mean         := (AverageHeight * MaxHeight) DIV 2;

 for Smpl:= 1 to NoOfSamples DO
  begin
   WndStart:= Smpl + 1;
   if WndStart > NoOfSamples then WndStart  := NoOfSamples;
   if WndStart < 1           then WndStart:= 1;

   WndEnd  := Smpl + 350;
   if WndEnd > NoOfSamples then WndEnd  := NoOfSamples;

   WndAverage:= AveragePeakHeightA(WndStart, WndEnd);
   if WndAverage = 0
   then WndAverage:= 1;

   if WndAverage < Mean
   then Chroma[Smpl].HeightA:= Chroma[Smpl].HeightA * RoundEx(Mean / WndAverage / 100)
   else EmptyDummy;
  end;
end;









function TCubeAbstractQV.AverageDomeHeight(CONST DomeStart, DomeEnd: Integer; Trace: RTrace): Cardinal;  { Media peak-urilor pt toata bucata de trace specificata }
VAR
   summ, Diff: Cardinal;
   i: Integer;
begin
 summ:= 0;
 for i:= DomeStart to DomeEnd
  DO summ:= summ+ Trace.Height[i];

 Diff:= DomeEnd-DomeStart;

 if Diff <= 0
 then Result:= 0
 else Result:= summ DIV Diff;
end;


function TCubeAbstractQV.MaxDomeHeight(CONST DomeStart, DomeEnd: Integer; Trace: RTrace): Cardinal;       { Highest peak in the specified dome }
VAR i: Integer;
begin
 Result:= 0;

 for i:= DomeStart to DomeEnd DO
  if Trace.Height[i] > Result
  then Result:= Trace.Height[i];
end;


function TCubeAbstractQV.MaxDomeHeight(CONST Sample: Integer; Trace: RTrace): Cardinal;
VAR DomeStart, DomeEnd: Integer;
begin
 DomeStart:= ComputeDomeStart (Sample, Trace);
 DomeEnd  := ComputeDomeEnd   (Sample, Trace);
 Result   := MaxDomeHeight(DomeStart, DomeEnd, Trace);
end;



{--------------------------------------------------------------------------------------------------
   QV - Average QV
--------------------------------------------------------------------------------------------------}
function TCubeAbstractQV.AverageQV: Integer;                                                       { Media QV-urilor pt bucata buna din secv }
VAR i: Integer;
begin
 Assert(FQVExist, 'This sample has no QV info.');

 Result:= 0;
 for i:= GoodQVStart to LastGoodBase
   DO Result:= Result+ CellQV[i];
 Result:= Result div NoOfGoodBases;
end;


function TCubeAbstractQV.AverageQVAll: Integer;                                                    { Media QV-urilor pt toata secv (low qv ends included) }
VAR i: Integer;
begin
 Assert(QVExist, 'This sample has no QV info.');

 Result:= 0;
 for i:= ctCellsIndex to NoOfBases
   DO Result:= Result+ CellQV[i];
 Result:= Result div NoOfBases;
end;





{--------------------------------------------------------------------------------------------------
   UTILS
--------------------------------------------------------------------------------------------------}
procedure TCubeAbstractQV.FreeTraces;                                                              { Free RAM }
begin
 SetLength(TraceA.Height, 0);
 SetLength(TraceC.Height, 0);
 SetLength(TraceG.Height, 0);
 SetLength(TraceT.Height, 0);
end;


procedure TCubeAbstractQV.ExtractTraces;                                                           { Extracting the 4 traces when I need them will be too slow. Instead I precalculate them and keep them in 4 local variables }
begin
 ExtractTrace('A', TraceA);
 ExtractTrace('C', TraceC);
 ExtractTrace('G', TraceG);
 ExtractTrace('T', TraceT);
end;


function TCubeAbstractQV.TraceNo2Trace(TraceNumber: Integer): RTrace;
begin
 case TraceNumber of
  1: Result:= TraceA;
  2: Result:= TraceC;
  3: Result:= TraceG;
  4: Result:= TraceT;
  else Result.Height:= NIL;                                                                        { Trace is NIL for non ACGT bases because it is difficult to tell exactly to which trace that base belongs. Theoretically N bases belong to all 4 traces, and UIPAC bases belong to one or more traces. But this will complicate the code toomuch }
 end;
end;

function TCubeAbstractQV.Base2Trace(Base: TBase): RTrace;
begin
  case Base of
   'A': Result:= TraceA;
   'C': Result:= TraceC;
   'G': Result:= TraceG;
   'T': Result:= TraceT;
   else Result.Height:= NIL;                                                                       { Trace is NIL for non ACGT bases because it is difficult to tell exactly to which trace that base belongs. Theoretically N bases belong to all 4 traces, and UIPAC bases belong to one or more traces. But this will complicate the code toomuch }
  end;
end;


procedure TCubeAbstractQV.MakeQVFromProbability;                                                   { The SCF does not have a single QV value. Instead it keeps 4 fields for each base (see TBase_v30). I consider that the highest value is the final QV }
VAR temp, CurBase: Integer;
begin
  for CurBase:= ctCellsIndex to NoOfBases DO                                                       { ambele matrici (Chroma si TempQV) sunt indexate in 1 }
   begin
    temp:= Find_Max(CellProbab[CurBase].A, CellProbab[CurBase].C, CellProbab[CurBase].G);
    CellQV[CurBase]:= Max(temp, CellProbab[CurBase].T);
   end;
end;





function TCubeAbstractQV.DomeThreshold: Real;                                                      { If a peak dome is under x% of AverageDomeArea then ignore it (consider it background) }
begin
 Result:= (SmallDomeThreshold / 100)* AverageDomeArea
end;





END.
