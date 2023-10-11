UNIT CubeBaseSNP;

{=============================================================================================================
 Gabriel Moraru
 2016.07
==============================================================================================================

 Object capabilities:
    + Double peak detection

 ----------------------------------------
  Parameters for SNP mutation detection
 ----------------------------------------
   ScanDistance: integer;                                                                          { Distance (between referece base and its neighbors) to scan for mutations. For 100% the program will scan for mutations a wider area. This will increase the chance to find mutations but it will also increase the chance to find false positives. For  10% the program will scan a narrow area. Default: 50%
   AreaRatio   : Double;                                                                           { Ratio (%) between polymorphic peak area and reference peak area. For example, a value of 10 means that the polymorphic peak area should be at least 10% of reference peak area. If the polymorphic peak area is under this value it is ignored. Default: 50%.
   OverlapArea : Integer;                                                                          { Overlap ratio (%) between polymorphic peak area and reference peak area. Ideally the peaks should fully overlap (the polymorphic peak should be under the reference peak). If the ratio is under this value, the mutation is ignored. Default: 70%.

  Test samples: d:\Biology\SAMPLES - Mutations\mutation between sequences\Staden mutations\

=============================================================================================================}

INTERFACE

USES System.SysUtils,
     ccCore, ccINIFile, CubeBase, CubicDNA, CubeBaseQv, clRamLog, ccRichLog;

TYPE
 TCubeAbstractSnp = class(TCubeAbstractQV)
  private
  protected
// procedure resetSnpOrig;
   function  DetectDoublePeak        (CONST BasePos: Integer; OUT OutTraceColor: TBase): Boolean;
   {DOME}
   function  ComputeDomeIntersection (CONST iFrom, iTo: Integer; MutTrace, MainTrace: RTrace): Integer;
   function  ComputeDomeStartEx      (CONST Sample : Integer; Trace: RTrace): Integer;                                         { Find where the 'dome' (with the center at Sample sample) starts/ends }  { Before I call this, I need to call ReasignPeaks to recompute dome center because ABI files are badly written }
   function  ComputeDomeEndEx        (CONST Sample : Integer; Trace: RTrace): Integer;
   function  ComputeDomeAreaEx       (CONST BasePos: Integer; Trace: RTrace; OUT DomeStart, DomeEnd: Integer): Integer;
  public
   EngSNP: REngSNP;
   TotalMutations: Integer;                                                                                                     { Total number of double peaks detected. Calculated by DetectDoublePeaks }
   constructor Create(aLog: TRamLog);
   procedure Clear; override;
   procedure DetectDoublePeaks;                                                                                                 { Called by TAsmJobAbst.LoadSample }  { Based on reference dome area/polymorphic dome area. See article: New strategy to detect single nucleotidepolymorphisms }
   function  GetPeakHeightRatio(const BasePos: Integer): string;
 end;


IMPLEMENTATION
USES System.Math, ccBinary, cmMath, Buckets;






{===============================================================================
   CREATE
===============================================================================}
constructor TCubeAbstractSnp.Create;                                                               { Unless you are careful, the object might not be fully constructed when the method is called. To avoid any problems, you should override the AfterConstruction method and use that for any code that needs to wait until the object is fully constructed. If you override AfterConstruction, be sure to call the inherited method, Top. }
begin
 inherited Create(aLog);
 SetDefaultSnpParam(EngSNP);                                                                       { DEFAULT MUTATION }
end;

procedure TCubeAbstractSnp.Clear;
begin
 inherited Clear;                                                                                  { Trebuie sa apelez neaparat inherited. Atunci cand apelez clear in cod, trebuei sa curat si parintii acestui obiect }
end;





{--------------------------------------------------------------------------------------------------
  DetectDoublePeaks

  Find double peaks based on reference dome area/polymorphic dome area (See article: New strategy to detect single nucleotidepolymorphisms)
  Applies only on the high quality (non-trimmed) region of the sequence.
  If a secondary peak (base) is detected then it is combined with the original base.
  The result will be an IUPAC base. It also sets the TotalMutations global var.
  Recommended to call RecallNBases before calling DetectDoublePeaks
  Called by TAsmJob.StartAssembly.
--------------------------------------------------------------------------------------------------}
procedure TCubeAbstractSnp.DetectDoublePeaks;
VAR
   CurBase, NewBase: TBase;
   iCurBase: Integer;
begin
 if NOT HasChroma
 then RAISE Exception.Create('Sample has no chromatogram: '+ ScreenName) at @TCubeAbstractSnp.DetectDoublePeaks;

 TotalMutations:= 0;

 { Reset bases to their original value }
 for iCurBase:= ctCellsIndex to NoOfBases DO
   if CellIsIUPAC(iCurBase)                                                                             { If the SnpOrig field is already #0 then I have no reason to change it }
   then Base[iCurBase]:= Cell(iCurBase).BaseOrig;                                                            { If I have a SNP in this point, I set it back to the original base }


 ExtractTraces;                                                                                         { Extracting the 4 traces when I need them will be too slow. Instead I precalculate them and keep them in 4 local variables }
 TRY
  { Peak identification for the 4 traces }
  for iCurBase:= GoodQVStart+1 to LastGoodBase-1 DO
   begin
    CurBase:= UpCase(Base[iCurBase]);

    { Skip N bases }
    if CharInSet(CurBase, DNA_GapN)
    then Continue;                                                                                      { If I still have N bases I just ignore them }

    { Do not recompute polimorphic (IUPAC) bases }
    if CharInSet(CurBase, Ambiguity) then
      begin
       Inc(TotalMutations);
       Continue;
      end;

    { Main algorithm }
    if DetectDoublePeak(iCurBase, NewBase) then
     begin
      Inc(TotalMutations);
      Base[iCurBase]:= MakeIUPAC(CurBase, NewBase);
     end;
   end;
 FINALLY
  FreeTraces;
 END;
end;


function TCubeAbstractSnp.DetectDoublePeak (CONST BasePos: Integer; OUT OutTraceColor: TBase): Boolean;   { Returns true if a base different than the current baser was found at this position }
VAR
   CurBase: TBase;
   MaxMutantArea: Integer;                                                                              { The highest area found between any of the 4 traces }
   LocalMutantArea: Integer;                                                                            { The highest area found between all mutants located under reference dome }
   RefTrace, MutantTrace: RTrace;
   RefBaseArea, CurTrace, ScanStart, ScanEnd, p, iOverlapArea: Integer;
   MainDomeStart, MainDomeEnd, MutantDomeStart, MutantDomeEnd: Integer;
   Mutant2RefArea, OverlapProc: Single;
begin
 Assert(BasePos>= ctCellsIndex);
 Assert(BasePos<= NoOfBases);

 Result:= FALSE;
 OutTraceColor:= noBase;
 CurBase:= UpCase(Base[BasePos]);

 Assert(CharInSet(CurBase, DNA_ACGT), 'Base at pos '+ IntToStr(BasePos)+ '[' + CurBase + ']+ should not be IUPAC!');

 RefTrace:= Base2Trace(CurBase);

 { Compute dome area for reference base }
 RefBaseArea:= ComputeDomeAreaEx(BasePos, RefTrace, MainDomeStart, MainDomeEnd);

 { Ignore tiny peaks }
 if RefBaseArea <= DomeThreshold                                                                   { If a peak dome is under x% of AverageDomeArea then ignore it (consider it background) }
 then EXIT(FALSE);                                                                                 { This is also 'Division by zero' protection }

 { Calculate range arround this base where to scan for mutations }
 ScanStart:= Base2Smpl[BasePos]- round( ProcentNormal(EngSNP.ScanDistance, DistanceBetween(BasePos- 1, BasePos))); // / ScanDistance);
 ScanEnd  := Base2Smpl[BasePos]+ round( ProcentNormal(EngSNP.ScanDistance, DistanceBetween(BasePos, BasePos+ 1)));

 { Peak identification for all three traces (self is autoexcluded) }
 MaxMutantArea:= 0;                                                                                { Find Max area among all 4 traces }
 for CurTrace:= 1 to 4 DO
  begin
   { Obtain current trace in a temp array to work on it }
   MutantTrace:= TraceNo2Trace(CurTrace);

   { Ignore self }
   if CurBase= MutantTrace.Color then Continue;                                                    { If I investigate peak A don't look for background peaks in its own trace }

   { For all samples in the selected trace }
   for p:= ScanStart+1 to ScanEnd-1 DO
    begin
       { Find peaks }
       if (MutantTrace.Height[p-1] <= MutantTrace.Height[p]) AND (MutantTrace.Height[p] >  MutantTrace.Height[p+1])
       OR (MutantTrace.Height[p-1] <  MutantTrace.Height[p]) AND (MutantTrace.Height[p] >= MutantTrace.Height[p+1]) then
        begin
          { A peak was found. Find its area }
          LocalMutantArea:= ComputeDomeArea(p, MutantTrace, MutantDomeStart, MutantDomeEnd);

          { Compute intersection: How much of this dome is under the main dome? }
          if (MutantDomeStart> MainDomeStart) AND (MutantDomeEnd< MainDomeEnd)                     { Mutant is totally included in the main dome }
          then iOverlapArea:= ComputeDomeIntersection(MutantDomeStart, MutantDomeEnd, MutantTrace, RefTrace)
          else
          if (MutantDomeEnd> MainDomeStart) AND (MutantDomeEnd<= MainDomeEnd)                      { Mutant is on the left side of the main dome }
          then iOverlapArea:= ComputeDomeIntersection(MainDomeStart, MutantDomeEnd, MutantTrace, RefTrace)
          else
          if (MutantDomeStart< MainDomeEnd) AND (MutantDomeStart>= MainDomeStart)                  { Mutant is on the right side of the main dome }
          then iOverlapArea:= ComputeDomeIntersection(MutantDomeStart, MainDomeEnd, MutantTrace, RefTrace)
          else
          if (MutantDomeStart<= MainDomeStart) AND (MutantDomeEnd>= MainDomeEnd)                   { Mutant is passign under the main dome but its ends are outside the main dome OR start exactly in the same position }
          then iOverlapArea:= ComputeDomeIntersection(MainDomeStart, MainDomeEnd, MutantTrace, RefTrace)
          else iOverlapArea:= 0;

          { Compare this peak area with the previous ones. Retain the the highest value among all mini domes found in the current trace }
          OverlapProc:= ProcentRepresent(iOverlapArea, LocalMutantArea);
          Mutant2RefArea:= (LocalMutantArea / RefBaseArea) * 100;

          if  (LocalMutantArea > MaxMutantArea)
          AND (OverlapProc     >= EngSNP.OverlapArea)                                              { and x% of this peak is under the main peak }       { If the ration between the reference peak area and mutant peak area is smaler than OverlapArea, than the mutant peak is ignored }
          AND (Mutant2RefArea  >= EngSNP.AreaRatio)
          then
           begin
            MaxMutantArea:= LocalMutantArea;
            OutTraceColor:= MutantTrace.Color;
            Result:= TRUE;                                                                         { Return the highest value among all 4 traces }

            { Store them for later use to show them in RamLog }
            CellsMX[BasePos].SnpAreaRatio:= IntToByte(RoundEx(Mutant2RefArea));
            CellsMX[BasePos].SnpOverlap  := IntToByte(RoundEx(OverlapProc));
           end;
        end;
    end;
  end;{END CurTRACE}
end;




function TCubeAbstractSnp.GetPeakHeightRatio(CONST BasePos: Integer): string;                      { Shown in RamLog }
VAR Sample: Integer;
    NewBase, OrigBase, IUPAC: TBase;
    Trace: RTrace;
begin
 Assert(CharInSet(Base[BasePos], ambiguity));
 Sample:= Base2Smpl[BasePos];

 { Split the IUPAC base at this position in its two components }
 OrigBase:= Cell(BasePos).BaseOrig;
 IUPAC   := Base[BasePos];

 NewBase := ExtractBase (OrigBase, IUPAC);                                                         { This is the new detected base }

 ExtractTrace(NewBase, Trace);                                                                     { Extracting the 4 traces when I need them is slow! }
 Result:= IntToStr(MaxDomeHeight(Sample, Trace));

 ExtractTrace(OrigBase, Trace);                                                                    { Extracting the 4 traces when I need them is slow! }
 Result:= Result+ ':'+ IntToStr(MaxDomeHeight(Sample, Trace));

 SetLength(Trace.Height, 0);
end;
















{--------------------------------------------------------------------------------------------------
   DOME START/END
--------------------------------------------------------------------------------
 Find where the 'dome' (with the center at Sample sample) starts/ends.
 Before I call this, I need to call ReasignPeaks to recompute dome center because ABI files are badly written.
 If Result= -1 it means that that the end of the dome cannot be found          }

CONST BackgroundNoise= 50;

function TCubeAbstractSnp.ComputeDomeStartEx(CONST Sample: Integer; Trace: RTrace): Integer;       { The 'ex' version computes also the plateau (it includes also 1/2 of the plateau between two bases in the dome) }
VAR
   Smpl, MidDistance: Integer;                                                                     { samples }
   PlateauHeight: Integer;
begin
 Result:= ComputeDomeStart(Sample, Trace);
 if Result= -1 then EXIT;

 { Identify the plateau region }
 PlateauHeight:= Trace.Height[Result];
 if PlateauHeight> BackgroundNoise then                                                            { plateau is considered only the region between two peaks (of the same color/trace) which height is <> than zero (BackgroundNoise) }
  for Smpl:= Result-1 downto ctChromaIndex DO
   if Trace.Height[Smpl] <> PlateauHeight then
    begin
     MidDistance:= (Result - (Smpl+1)) DIV 2;                                                      { Find plateau midpoint } {DIV = the result is rounded in the direction of zero to the nearest Integer }
     Result:= Result - MidDistance;                                                                { The plateau midpoint is considered as the region limit }
     Break;
    end;
end;


function TCubeAbstractSnp.ComputeDomeEndEx(CONST Sample: Integer; Trace: RTrace): Integer;
VAR
   Smpl, MidDistance: Integer;
   PlateauHeight: Integer;
begin
 Result:= ComputeDomeEnd(Sample, Trace);
 if Result= -1 then EXIT;

 PlateauHeight:= Trace.Height[Result];
 if PlateauHeight> BackgroundNoise then                                                            { plateau is considered only the region between two peaks (of the same color/trace) which height is <> than zero (BackgroundNoise) }
  for Smpl:= Result+1 to NoOfSamples DO
   if Trace.Height[Smpl] <> PlateauHeight then
    begin
     MidDistance:= ((Smpl-1) - Result) DIV 2;
     Result:= Result + MidDistance;
     Break;
    end;
end;








{--------------------------------------------------------------------------------------------------
   DOME AREA
--------------------------------------------------------------------------------------------------}
function TCubeAbstractSnp.ComputeDomeAreaEx(CONST BasePos: Integer; Trace: RTrace; OUT DomeStart, DomeEnd: Integer): Integer;  { 'Trace' is the trace to which this base belongs }
VAR
   Sample, Smpl: Integer;
begin
 { Find margins }
 Sample:= Base2Smpl[BasePos];
 DomeStart := ComputeDomeStartEx (Sample, Trace);
 if (DomeStart<= 0) then EXIT(0);

 DomeEnd:= ComputeDomeEndEx   (Sample, Trace);
 if (DomeEnd<= 0) then EXIT(0);

 Assert(DomeStart <= DomeEnd);

 { Compute area }
 Result:= 0;
 for Smpl:= DomeStart to DomeEnd DO
  Result:= Result+ Trace.Height[Smpl];
end;




{--------------------------------------------------------------------------------------------------
   DOME INTERSECT/REUNION
--------------------------------------------------------------------------------------------------}
function TCubeAbstractSnp.ComputeDomeIntersection(CONST iFrom, iTo: Integer; MutTrace, MainTrace: RTrace): Integer;
VAR Smpl, MinIntersect: Integer;
begin
 Result:= 0;
 Assert(iFrom<= iTo);
 Assert(iFrom>= ctChromaIndex);
 Assert(iFrom<= NoOfSamples);

 for Smpl:= iFrom to iTo DO
  begin
   MinIntersect:= Min(MutTrace.Height[Smpl], MainTrace.Height[Smpl]);
   Result:= Result+ MinIntersect;
  end;
end;





END.
