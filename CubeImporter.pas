
UNIT CubeImporter;

{===============================================================================
 Heracle BioSoft SRL
 2016.10.12

 Object capabilities:
   + Import from  SCF, FASTA, SEQ, ABI
   + Export to    SCF, FASTA, SEQ
===============================================================================}

INTERFACE

USES
   Winapi.Windows, System.SysUtils, System.Math, Vcl.Graphics,
   CubicDNA, ScfRead, ScfBase, ReadAbi, ReadFasta, ParseGenBankR, CubeBase, CubeBaseSNP, ccCore, ccINIFile, clRamLog, ccRichLog;

TYPE
 TCubeImport = class(TCubeAbstractSnp)
  private
   FParentVer : string;
   FUseIntQvAlgo: Boolean;                                                                                    { if true, I will always use my own algorithm to calculate/recalculate QVs. If false, I will use the QV calculated by Phred. If no QV data exista I use internal algorithm to calculate the QV }
   FRecallNPeaks: Boolean;
  protected
   function Boot: Boolean;                                                                                    { Post processing after laoding the object. This function also calls TrimmEngine }
   function cloneAndTrim_: TCubeImport; //Check: DO I REALLY NEED IT? Get rid of it                           { Return a copy of the curent cube with low quality ends (determined with TrimEngine) trimmed }
  public
   CommentOrig: string;                                                                                       { original comment. Used only for SCF and ABI }
   OrigLength : Integer;                                                                                      { Original seq length, before trimming }
   ChromaDisplay: TObject;                                                                                    {ASSOCIATED DISPLAY}        { pointer to a associated TChromaDisplay object }

   constructor Create(aLog: TRamLog);
   procedure   Clear;   override;

   {LOG & MESSAGES}
   function  PropertiesWithName  : string;                                                                    { before I call this I have to apply TrimEngine }
   function  Properties          : string;                                                                    { before I call this I have to apply TrimEngine }
   function  ShowEstimatedQuality: string;
   function  ShowBases           : string;
   function  ShowGoodBases       : string;
   function  ShowBasesAfterTrim  : string;                                                                    { spune cat % din baze au ramas dupa TrimEngine }
   function  ShowQVExist         : string;
   function  TrustedBasesPercent: Real;                                                                       { How many trusted bases (bases with QV over GoodQVTresh) are there? }

   {IMPORT}
   function  Import(CONST FileName: string): Boolean;                                                         { Fasta and GBK files are not supported because they may contain more than one sample! }
   function  AssignGBK     (CONST Gbk : TGbkObj): Boolean;
   function  AssignScf     (CONST SCF : TScfObj): boolean;
   function  AssignAbi     (CONST ABI : TAbiObj; EditedField: boolean= True): Boolean;                        { parametrul EditedFiled= arata de unde sa extraga informatia: din campul original sau din campul EDITED }
   function  AssignSecv    (CONST SECV: TFastaObj): Boolean;
   property  Version: string  read FParentVer  Write FParentVer;                                              { ce versiune a avut obiectul parinte. Exemplu, pt SCF: '3.00' }

   {EXPORT}
   function  ExportBack    (CONST Clean: Boolean; OUT BasesLeft: Integer): Boolean;                           { Save the object back to disk, in its original format. This will overwrite the parent sample }
   function  ExportAs      (CONST aFileName: string; Clean: Boolean; OUT BasesLeft: Integer): Boolean;        { Autodetecteaza tipul secventei din numele fisierului. QvAware=True transfera numai cu bazele QV bun }
   procedure ExportAsFASTA (CONST                    Clean: Boolean; OUT BasesLeft: Integer); overload;       { uses cube's name and changes its extension to 'fasta' }
   procedure ExportAsFASTA (CONST FileName : string; Clean: Boolean; OUT BasesLeft: Integer); overload;
   procedure ExportAsSecv  (CONST NewName  : string; Clean: Boolean; OUT BasesLeft: Integer);
   procedure ExportToSCF   (CONST aFileName: string; Clean: Boolean);            overload;                    { This will remove all GAPs }
   function  ExportToSCF   (CONST                    Clean: Boolean): string;    overload;                    { Returns the name used to save the file }
   function  AsFasta       (CONST aFileName: string; Clean: Boolean): TFastaObj; overload;                    { Return a TFastaObj object built from this cube }
   function  AsFasta       (CONST                    Clean: Boolean): TFastaObj; overload;                    { Same as above but the user doesnt have to provide a file name. It is automatically assumed from Cube's name }

   {BUILD}
   procedure BuildFromScratch (CONST sBases: BaseString; CONST sFullName, sComments: string);                 { Build a cube from the given parameters }
   property  RecomputeQVs  : Boolean read  FUseIntQvAlgo write FUseIntQvAlgo default TRUE;                    { if true, I will always use my own algorithm to calculate/recalculate QVs. If false, I will use the QV calculated by Phred. If no QV data exista I use internal algorithm to calculate the QV }
   property  RecallNPeaks  : Boolean read  FRecallNPeaks write FRecallNPeaks default TRUE;

   {DEBUG}
   function  DebugToBmp    (ShowLinks: Boolean; YDivider, XMultiplier: integer): TBitmap;                     { creates a bitmap that shows the internal laying of the chromatogram }
   procedure DebugToBmpSave(SaveTo: string; ShowLinks: Boolean; CONST YDivider: Integer= 10; CONST XMultiplier: Integer= 5);   { creates a bitmap that shows the internal laying of the chromatogram. 10 is a good default value for YDivider. if Path is empty, the routine will autogenerate the path }
   procedure DebugToTxt    (Path: string);
   procedure FakeChroma    (PeakHeight: integer; Randomiz: boolean);                                          { FAKE A CHROMA - PeakHeight should be about set at 700. ATENTIE: trebuie sa apelez "Randomize" inainte a a apela "FakeChroma" }
   procedure SaveChromaAsBmp(CONST Nume: string);
 end;


IMPLEMENTATION
USES
  ccAppData, ccIO, cmMath, SCFWrite;





{===============================================================================
   CONSTRUCTOR / DESTRUCTOR
===============================================================================}

constructor TCubeImport.Create;
begin
 inherited Create(aLog);
 FUseIntQvAlgo := TRUE;                                                                            { If true, it means that the internal algorithm was used to compute the QVs }
 FRecallNPeaks := TRUE;
end;



{ Call Clear when:
     * create an abstract cube (TCubeBase, TCubeImport, TCubObj), then call Clear immediately after Create.
     * create a persistent cube (TCubObjEx), then call Clear when I call LoadFromFile or Assign (if I ever implement Assign). }
procedure TCubeImport.Clear;
begin
 inherited Clear;                                                                                  { I have to call inherited. When I call clear in the code, I must also clear the parents of this object }
 Version    := '';                                                                                 { what version the parent object had. Example, for SCF: '3.00' }
 CommentOrig:= '';
 OrigLength := 0;
end;




function TCubeImport.Boot: Boolean;                                                                { Post processing after laoding the object. This function also calls TrimmEngine }

  procedure CheckQVExists;
  VAR i: Integer;
  begin
   if NOT HasChroma
   then FQVExist:= FALSE
   else
     for i:= 1 TO NoOfBases DO
      if CellQV[i] > 1 then                                                                        { I chose 1 instead of 0, to treat case where QV field was:  0, 1, 0, 1...  ) }
       begin
        FQVExist:= TRUE;
        Break;
       end;                                                                                        { ReturnInfo_:= 'File name: '+ ShortName +CRLF+ 'File path: '+ ExtractFilePath(ParentFileName) +LBRK+ ABI.ReturnInfo1}
  end;

begin
 Result:= TRUE;

 { CHECK AGAINST TINY CHROMAS }                                                                    { see 'amoA pCC1 094.scf' }
 if NoOfBases<= ctMinimimSequence then
  begin
    RamLog.AddWarn('This sample has less than 10 bases!');
    EXIT(FALSE);
  end;

 { QV exists? }
 CheckQVExists;

 { BASES }
 DirtyBases:= TRUE;
 OrigLength:= NoOfBases;

 if HasChroma AND (NoOfBases> 0) then
  begin                                                                                            { Make sure that all samples start and end at the ground }
   Chroma[ctChromaIndex].HeightA:= 1;
   Chroma[ctChromaIndex].HeightC:= 1;
   Chroma[ctChromaIndex].HeightG:= 1;
   Chroma[ctChromaIndex].HeightT:= 1;
   Chroma[NoOfSamples]  .HeightA:= 1;
   Chroma[NoOfSamples]  .HeightC:= 1;
   Chroma[NoOfSamples]  .HeightG:= 1;
   Chroma[NoOfSamples]  .HeightT:= 1;
  end;

 ReasignPeaks;                                                                                     { For some strange reasons, in ABI files, the pointers between the base and the peak are placed few points BEFORE the top of the dome (the real peak). So I have to recalculate their position }

 { CALCULATE QVs }
 if HasChroma then
  begin
   if RecallNPeaks
   then RecallNBases;                                                                              { Recall the N bases }

   if NOT QVExist
   then CalculateQVs                                                                               { No QV included. Force CalculateQVs }
   else if RecomputeQVs                                                                            { If 'Prefer internal base caller' then force CalculateQVs }
        then CalculateQVs;
  end;

 buildBases;

 { TRIM ENGINE }
 TrimEnds;
end;







{===============================================================================
   IMPORTER

   This calls TrimEnds and DetectVectors
   Note: Fasta and GBK files might contain more than one sample! Only the first sample will be loaded
===============================================================================}
function TCubeImport.Import(CONST FileName: string): Boolean;
VAR
  unSCF: TScfRead;
  unABI: TAbiObj;
  Fasta: TFastaObj;
begin
 Result:= FALSE;                                                                                   { Must be initialized otherwise it keeps in STACK the value on the previous run }

 { LOAD SCF }
 if IsScf(FileName) then
  begin
   unSCF:= TScfRead.Create(RamLog);
   TRY
     if unSCF.LoadFromFile(FileName)
     then Result:= AssignScf(unSCF);
   FINALLY
     FreeAndNil(unSCF);
   END
  end
 else

 { LOAD ABI }
 if IsAbi(FileName) then
  BEGIN
   unABI:= TAbiObj.Create(RamLog);
   TRY
     if unABI.LoadFromFile(FileName, FALSE)
     then Result:= AssignAbi(unABI, TRUE);
   FINALLY
     FreeAndNil(unABI);
   END
  END
 else

 { LOAD FASTA/SEQ/TXT/GBK }
 if IsPlainText(FileName) then
  BEGIN
   Fasta:= TFastaObj.Create(RamLog);
   TRY
     if Fasta.LoadFromFile(FileName, FALSE)
     then Result:= AssignSecv(Fasta);
   FINALLY
     FreeAndNil(Fasta);
   END
  END
 else
    RamLog.addError('Unsupported file: '+ FileName);                                             { Only for those with chromatograms }

 //del - this is done in Boot         -> buildBases;
end;






{===============================================================================
   ASSIGN ABI
===============================================================================}

function TCubeImport.AssignAbi(CONST ABI: TAbiObj; EditedField: boolean= True): Boolean;  { the EditedFiled= parameter shows where to extract the information: from the original field or from the EDITED field }
VAR
   ProbRec: RBaseProb;

  procedure ImportBases(Side: TSide);                                                 { Import Bases and QVs }
  VAR
     CurBase: Integer;
     Ptr2Smpl: Integer;
     QV: Byte;
  begin
    NoOfBases := Side.NrOfBases;                                                      { CellMX is indexed in 1 not 0 }
    CellMxReset;                                                                      { I must call this AFTER I set NoOfBases }

    Assert(Length(Chroma)= NoOfSamples+1, 'Invalid chroma size.');

    for CurBase:= ctCellsIndex to NoOfBases DO
     begin
      { Import BASE }
      Base[CurBase]:= tbase(Side.BaseArray[CurBase]);                                 { set Base, but it doesn't call 'Changed' }

      { Import POINTERS }
      Ptr2Smpl:= Side.Ptr2Smpl[CurBase];                                              { must be equal to the number of bases, so the association is from BASES TO SAMPLES }

      if Ptr2Smpl >= Length(Chroma)                                                   { This fixes the Russell case}
      then Ptr2Smpl:= Length(Chroma)-1;

      Chroma[Ptr2Smpl].Ptr2Base := CurBase;                                           { BASE POS }
      Base2Smpl[CurBase]:= Ptr2Smpl;
     end;

    { Import QV }
    if Length(Side.QvMX) > 0 then                                                     { Some ABI files will not have QV info, so the QvMx will be empty }
     begin
      Assert(NoOfBases= Length(Side.QvMX), 'NoOfBases different than QVCount!');
      for CurBase:= ctCellsIndex to NoOfBases DO
       begin
        QV:= Side.QvMX[CurBase-IndexDiff];
        CellQV[CurBase]:= QV;

        { Convert ABI QV to SCF probabilities }                                       { SCF format has this feature but ABI format don't }
       if UpCase(Base[CurBase])= 'N' then                                             { For N, all bases have same probability (same QV) }
        begin
         ProbRec.A:= QV;
         ProbRec.C:= QV;
         ProbRec.G:= QV;
         ProbRec.T:= QV;
        end
       else
        begin
         FillChar(ProbRec, SizeOf(ProbRec), #0);                                       { Necessary because records are not automatically initialized! }
         case UpCase(Base[CurBase]) of
          'A': ProbRec.A:= QV;
          'C': ProbRec.C:= QV;
          'G': ProbRec.G:= QV;
          'T': ProbRec.T:= QV;
         end;
        end;

        CellsMX[CurBase].CellProb:= ProbRec;
       end;
     end;
  end;


  { Import traces (into samples) }
  procedure ImportTraces(Trace: TAbiTrace);
  VAR j: Integer;
  begin
   case Trace.Name of
    'A': for j:= 0 to NoOfSamples-1 DO Chroma[j+1].HeightA:= Trace.points[j];
    'C': for j:= 0 to NoOfSamples-1 DO Chroma[j+1].HeightC:= Trace.points[j];
    'G': for j:= 0 to NoOfSamples-1 DO Chroma[j+1].HeightG:= Trace.points[j];
    'T': for j:= 0 to NoOfSamples-1 DO Chroma[j+1].HeightT:= Trace.points[j];
   end;
  end;


begin
 Clear;                                                                               { THIS IS MANDATORY in order to clear the previous data (for example pointerMX ) loaded in Cube. I tested the program without it and didn't worked. }
 TRY
  FileName:= ABI.FileName;
  Comment := '> '+ ScreenName;                                                        { ABI has no comment field so 'Comment' will be empty }     { Only FASTA files have relevant comments. For SCF and ABI we ignore the comments. Cristina}
  CommentOrig:= string(ABI.H.AbiTags);
  NoOfSamples:= ABI.NoOfSamples;                                                      { This MUST be ABOVE 'import bases' }

  { Import Bases and QVs }
  if EditedField
  then ImportBases(ABI.Side1)
  else ImportBases(ABI.Side2);

  { Import traces (into samples) }
  ImportTraces(ABI.Trace1);
  ImportTraces(ABI.Trace2);
  ImportTraces(ABI.Trace3);
  ImportTraces(ABI.Trace4);

  Result:= Boot;
 EXCEPT
   on E : Exception DO
    begin
     RamLog.AddError(E.Message + CRLF+ 'Invalid ABI file. Please report this error to us.');
     Result:= FALSE;
    end;
 END;
end;



{===============================================================================
   ASSIGN SCF
===============================================================================}

function TCubeImport.AssignScf(const SCF: TScfObj): Boolean;
VAR i, CurBase: Integer;
    ProbRec: RBaseProb;
    iPointer: Cardinal;
begin
 Clear;                                                                                            { THIS IS MANDATORY in order to clear the previous data (for example pointerMX ) loaded in Cube. I tested the program without it and didn't worked. }

 TRY
  ParentType    := bfSCF;
  NoOfSamples   := SCF.NrOfSamples;                                                                { reserve memory }
  NoOfBases     := SCF.NoOfBases;                                                                  { reserve memory }
  Version       := String(SCF.H.Version);
  FileName      := SCF.FileName;

  Comment       := '> '+ ScreenName;
  CommentOrig   := SCF.Comments;                                                                   { original comment. Used only for SCF and ABI }                                                            { Some scf files have no comment, others have something like this: NAME=CMO_59 LANE=1  SIGN=A=652,C=962,G=869,T=873 }{ Only FASTA files have relevant comments. For SCF and ABI we ignore the comments. Cristina}

  CellMxReset;

  for CurBase:= ctCellsIndex to NoOfBases DO                                                       { both matrices (Chroma and TempQV) are indexed in 1 }
   begin
    { Import BASE }
    Base[CurBase]:= tbase(SCF.BaseArray[CurBase].Base);

    { Import pointers }
    iPointer:= SCF.BaseArray[CurBase].Ptr2Smpl;
    if iPointer> SCF.NrOfSamples                                                                   { 'Amoa sample' special case }
    then iPointer:= 1;
    Chroma[iPointer].Ptr2Base := CurBase;
    Base2Smpl[CurBase]:= iPointer;

    { Import probabilities }
    ProbRec.A:= SCF.BaseArray[CurBase].prob_A;
    ProbRec.C:= SCF.BaseArray[CurBase].prob_C;
    ProbRec.G:= SCF.BaseArray[CurBase].prob_G;
    ProbRec.T:= SCF.BaseArray[CurBase].prob_T;

    CellsMX[CurBase].CellProb:= ProbRec;
   end;

  { Compute QVs for all bases }
  MakeQVFromProbability;

  { Import trace-urile }                                                    
  for i:= 0 to NoOfSamples-1 DO                                                                    { SCF is indexed in 0 }
   begin
    Chroma[i+ ctChromaIndex].HeightA:= SCF.TraceA[i];                                              { +1 because my Chroma matrix is indexed in 1 }
    Chroma[i+ ctChromaIndex].HeightC:= SCF.TraceC[i];
    Chroma[i+ ctChromaIndex].HeightG:= SCF.TraceG[i];
    Chroma[i+ ctChromaIndex].HeightT:= SCF.TraceT[i];
   end;

  { In v3.10 QV-urile sunt stocate ca Log10 }
  if Version>= '3.10'then
   for i:= 0 to NoOfSamples-1 DO
    begin
     CellQV [i+1] := round(-10 * Log10(CellQV[i+1]));
     if CellQV [i+1] < 1
     then CellQV [i+1] := 0;
    end;

  Result:= Boot;

 EXCEPT
   RamLog.AddError('SCF file cannot be imported. Please report this error to us.');
   Result:= FALSE;
 END;
end;







function TCubeImport.ExportToSCF(CONST Clean: Boolean): string;                                  { Returns the name used to save the file. }
begin
 Result:= ChangeFileExt(FileName, '.scf');
 ExportToSCF(Result, Clean);                                                                     { This will remove all GAPs }
end;


procedure TCubeImport.ExportToSCF(CONST aFileName: string; Clean: Boolean);                      { This will remove all GAPs }
VAR UnScf: TScfWriter;
    CurSmpl, CubBase, ScfBase, Ptr2Smpl: Integer;
    Duplicate: TCubeImport;
    Base: TBase;
begin
 if Clean AND FQVExist
 then Duplicate:= cloneAndTrim_                                                                  { CloneAndTrim returns a copy of the curent cube with low quality ends (determined with TrimEnds) trimmed. The caller has to free the cube. }
 else Duplicate:= Self;

 UnScf:= TScfWriter.Create(RamLog);
 TRY
  UnScf.NrOfSamples  := Duplicate.NoOfSamples;                                                   { This will make:  SetLength }
  UnScf.NoOfBases    := Duplicate.NoOfBases;                                                     { This will make:  SetLength(TempQV, H.NoOfBases+ IndexedIn1) }
  UnScf.H.Version    := '3.00';
  UnScf.Comments     := Duplicate.CommentOrig;
  UnScf.H.SampleSize := 2;
  UnScf.H.ClipLeft   := 0;                                                                       { Obsolete. We don't care. I just initialize them to 0 so when I write this field to disk it won't write a random value }
  UnScf.H.ClipRight  := 0;

  { Transfera bazele }
  ScfBase:= ctCellsIndex;
  for CubBase:= ctCellsIndex to Duplicate.NoOfBases DO
   begin
    Base:= Duplicate.Base[CubBase];

    if Base<> GAP
    then                                                                                         { Skip GAPs }
     begin
      Ptr2Smpl:= Duplicate.Base2Smpl[CubBase];
      Assert(Ptr2Smpl<= Duplicate.NoOfSamples, ' Ptr2Smpl= '+i2s(Ptr2Smpl));

      UnSCF.BaseArray[ScfBase].Ptr2Smpl:= Ptr2Smpl;
      UnSCF.BaseArray[ScfBase].Base:= Ansichar(Base);

      { Transfer SCF probabilities }                        { SCF supports per-base probabilities }
      if ParentType= bfSCF then
       begin
        UnSCF.BaseArray[ScfBase].prob_A:= Duplicate.CellProbab[CubBase].A;
        UnSCF.BaseArray[ScfBase].prob_C:= Duplicate.CellProbab[CubBase].C;
        UnSCF.BaseArray[ScfBase].prob_G:= Duplicate.CellProbab[CubBase].G;
        UnSCF.BaseArray[ScfBase].prob_T:= Duplicate.CellProbab[CubBase].T;
       end;

      { Transfer QVs }               { However, if a base is edited in DNA Baser editor (or by the internal base caller), the probability is stored in the QV field not in the Probabilities field. So, I need to stansfer also the QVs }
      case Duplicate.Base[CubBase] of
       'A': UnSCF.BaseArray[ScfBase].prob_A:= Duplicate.CellQV[CubBase];
       'C': UnSCF.BaseArray[ScfBase].prob_C:= Duplicate.CellQV[CubBase];
       'G': UnSCF.BaseArray[ScfBase].prob_G:= Duplicate.CellQV[CubBase];
       'T': UnSCF.BaseArray[ScfBase].prob_T:= Duplicate.CellQV[CubBase];
       else
         begin
           UnSCF.BaseArray[ScfBase].prob_A:= Duplicate.CellQV[CubBase];
           UnSCF.BaseArray[ScfBase].prob_C:= Duplicate.CellQV[CubBase];
           UnSCF.BaseArray[ScfBase].prob_G:= Duplicate.CellQV[CubBase];
           UnSCF.BaseArray[ScfBase].prob_T:= Duplicate.CellQV[CubBase];
         end;
      end;

      Inc(ScfBase);
     end
    else EmptyDummy;    { Skip GAPs }
   end;
  UnScf.NoOfBases:= ScfBase-1;                                                                 { Take into account the removed GAPs }

  { Transfera trace-urile }
  for CurSmpl:= 0 to Duplicate.NoOfSamples-1 DO                                                { In SCF ChromaMX e indexata in 0 }
   begin
    UnSCF.TraceA[CurSmpl]:= Duplicate.Chroma[CurSmpl+ ctChromaIndex].HeightA;
    UnSCF.TraceC[CurSmpl]:= Duplicate.Chroma[CurSmpl+ ctChromaIndex].HeightC;
    UnSCF.TraceG[CurSmpl]:= Duplicate.Chroma[CurSmpl+ ctChromaIndex].HeightG;
    UnSCF.TraceT[CurSmpl]:= Duplicate.Chroma[CurSmpl+ ctChromaIndex].HeightT;
   end;

  UnScf.SaveToFile(aFileName);

 FINALLY
  FreeAndNil(UnScf);
  if Duplicate<> Self          //del Clean AND FQVExist
  then FreeAndNil(Duplicate);
 END;

end;


function TCubeImport.cloneAndTrim_: TCubeImport;                                               { Returns a copy of the curent cube with low quality ends (determined with TrimEnds) trimmed. The caller has to free the cube. }
VAR
   b1, b2, SampleDiff, BaseDiff: Integer;
   FromBase, ToBase: Integer;
   FromSample, ToSample: Integer;
   Sample1, Sample2: Integer;
begin
 Assert(FQVExist);
 Result:= TCubeImport.Create(RamLog);

 { FROM/TO }
 CleanFromTo(FromBase, ToBase);                                                                { Shows from which to which base I have to cut in order to remove the vectors and the low quality ends }

 { Too short? }                                                                                { If the sequence left after trimming is to short then I force it to have at least 5 good bases. }
 if ToBase-FromBase < ctMinimimSequence then
  begin
   FromBase:= NoOfBases DIV 2;
   ToBase:= FromBase+ ctMinimimSequence;
   if ToBase > NoOfBases
   then ToBase:= NoOfBases; // RAISE Exception.Create('Chromatogram is too short! '+ FileName);
  end;

 FromSample:= Base2Smpl[FromBase] - 8;                                                         { 6 is half of average peak distance (12) but we take a bit more because some somes are larger than 12 samples. We need it in order to copy the left side of the dome. Otherwise we start copy samples wight in the middle of the dome }
 ToSample  := Base2Smpl[ToBase]   + 8;

 if FromSample < ctChromaIndex
 then FromSample:= ctChromaIndex;
 if ToSample> NoOfSamples
 then ToSample:= NoOfSamples;

 {cube PROPERTIES}
 Result.NoOfBases   := 1+ ToBase- FromBase;
 Result.NoOfSamples := NoOfSamples;                                                       { DUMMY CALL. I need it because in TCubeAbstract.GetPtr2Smpl I have an Assertion that blows. I make the REAL assignment few lines down down }
 Result.FQVExist    := QVExist;
 Result.Comment     := Comment;
 Result.Reversed    := Reversed;

 {BASES/QV}
 Result.EngTrim1    := EngTrim1;                                                          { the left end - for the moment it is used for both ends }
 Result.EngTrim2    := EngTrim2;                                                          { the right end }

 {CUT IT}
 b1 := 0;
 SampleDiff:= FromSample-1;
 BaseDiff  := FromBase-1;

 { Copy from specified bases }
 for b2:= FromBase to ToBase DO
  begin
    Inc(b1);
    Result.CellSet(b1, Cell(b2));
    { Refac pointerii }
    Result.Base2Smpl[b1]:= Result.Base2Smpl[b1]- SampleDiff;
  end;

 Sample1:= 0;
 for Sample2:= FromSample to ToSample DO
  begin
    Inc(Sample1);
    Result.Chroma[Sample1]:= Chroma[Sample2];
    { Refac pointerii }
    if Result.Chroma[Sample1].Ptr2Base > 0
    then Result.Chroma[Sample1].Ptr2Base:= Result.Chroma[Sample1].Ptr2Base- BaseDiff;
  end;

 DirtyBases:= TRUE;
 DetectVectors;
 Result.NoOfSamples:= 1+ ToSample- FromSample;

 Result.GoodQVStart:= 1;
 //Result.GoodQVEndS:= NoOfSamples;

 Result.TrimEnds;    { Calls GetGoodBases, GoodQVStart, GoodQVEnd, GoodQVStaSmpl, GoodQVEndSmpl si CHANGED }
end;










{ TRANSFER from plain sequence (SECV) }
function TCubeImport.AssignSecv(CONST SECV: TFastaObj): Boolean;
begin
 Result:= FALSE;
 Clear;                                                                                            { THIS IS MANDATORY in order to clear the previous data (for example pointerMX ) loaded in Cube. I tested the program without it and didn't worked. }

 if (SECV.NoOfBases<= 1) then
  begin
    RamLog.AddError('Invalid file: the sample is empty.');
    EXIT(FALSE);
  end;

 TRY
  { FILE }
  FileName  := SECV.FileName;                                                                      { Note: I should use SECV.ParentName here but in this case it won't see the names assigned to multi-samples by the GenerateVirtualNames (for multi-fasta files) }
  Comment   := SECV.Comment;
  IsPart    := SECV.IsPart;

  { cube }
  NoOfBases     := SECV.NoOfBases;
  CellMxReset;                                                                                     { Set pointers and QVs to 0 }
  NoOfSamples   := 0;
  FQVExist      := FALSE;

  { Importa bazele }
  setBases(SECV.BASES);

  Result:= Boot;
 EXCEPT
  on E: Exception
   DO RamLog.AddError('Error while trying to read the sample. '+ E.Message);
 END;
end;



{ Import from GBK }
function TCubeImport.AssignGBK(CONST Gbk: TGbkObj): Boolean;
begin
 Clear;                                                                                            { THIS IS MANDATORY in order to clear the previous data (for example pointerMX ) loaded in Cube. I tested the program without it and didn't worked. }

 if (Gbk.NoOfBases<= 1) then
  begin
    RamLog.AddError('Invalid file: the sample is empty.');
    EXIT(FALSE);
  end;

 TRY
  {FILE}
  FileName:= Gbk.FileName;
  Comment := Gbk.Comment;
  IsPart  := Gbk.IsPart;

  {cube}
  NoOfBases:= Gbk.NoOfBases;                                                                  { reserve memory }
  CellMxReset;

  NoOfSamples:= 0;
  FQVExist   := FALSE;

  {Import bases}
  setBases(Gbk.BASES);

  Result:= Boot;
 EXCEPT
   RamLog.AddError('Invalid GBK file. Please report this error to us.');
   Result:= FALSE;
 END;
end;






{===============================================================================
   EXPORT IT IN THE ORIGINAL FORMAT 
===============================================================================}
function TCubeImport.ExportBack(CONST Clean: Boolean; OUT BasesLeft: Integer): Boolean;            { Save the object back to disk, in its original format. This will overwrite the parent sample. Useful after the user changed some bases in it }
begin
 Result:= TRUE;
 case ParentType of
   bfSCF: Result:= ExportToSCF  (Clean)= '';                                                       { This will remove all GAPs }
   bfTXT: ExportAsSecv (FileName, Clean, BasesLeft);
   bfSEQ: ExportAsSecv (FileName, Clean, BasesLeft);
   bfFAS: begin Result:= TRUE; ExportAsFASTA(Clean, BasesLeft); end;
 else    
   Result:= FALSE; { bfABI, bfGBK, bfNone }
 end;
end;





{===============================================================================
   EXPORT TO A DIFFERENT FORMAT
===============================================================================}
function TCubeImport.ExportAs(CONST aFileName: string; Clean: Boolean; OUT BasesLeft: Integer): Boolean;  { Autodetecteaza tipul secventei din numele fisierului. QvAware=True transfera numai cu bazele QV bun }
begin
 Result:= TRUE;
 if IsScf(aFileName)
 then ExportToSCF(aFileName, Clean)                                                                 { This will remove all GAPs }
 else
   if IsPlainText(aFileName)
   then ExportAsSecv(aFileName, Clean, BasesLeft)
   else Result:= FALSE;                                                                            { nu am putut sa salvez pt ca nu recunosc extensia asta }
end;


procedure TCubeImport.ExportAsSecv(CONST NewName: string; Clean: Boolean; OUT BasesLeft: Integer); { 'NewName' is full path. Clean means "Remove low quality ends and vectors" }
VAR Secv: TFastaObj;
begin
 Secv:= AsFasta(NewName, Clean);
 TRY
  BasesLeft:= Secv.NoOfBases;
  Secv.Save_AutoDetect;
 FINALLY
  FreeAndNil(Secv);
 END;
end;






procedure TCubeImport.ExportAsFASTA(CONST Clean: Boolean; OUT BasesLeft: Integer);         { uses cube's name and changes its extension to 'fasta' }
begin
 ExportAsFASTA(ChangeFileExt(FileName, '.FASTA'), Clean, BasesLeft);
end;


procedure TCubeImport.ExportAsFASTA(CONST FileName: string; Clean: Boolean; OUT BasesLeft: Integer);
VAR Secv: TFastaObj;
begin
 Assert(SameText(ExtractFileExt(FileName), '.fasta'), 'File ext MUST be Fasta');

 Secv:= AsFasta(FileName, Clean);
 TRY
  BasesLeft:= Secv.NoOfBases;
  Secv.Save(TRUE);
 FINALLY
  FreeAndNil(Secv);
 END;
end;





function TCubeImport.AsFasta(CONST aFileName: string; Clean: Boolean): TFastaObj;                  { Return a TFastaObj object built from this cube. The caller needs to free it }
begin
 Result:= TFastaObj.Create(RamLog);

 { Remove low quality ends and vectors }   { It also removes GAPS }
 if Clean
 then Result.BASES := GoodBasesNoVector
 else Result.BASES := Bases;

 Result.Comment := ValidateComment(Comment);
 Result.FileName:= aFileName;
end;



function TCubeImport.AsFasta(CONST Clean: Boolean): TFastaObj;                                    { Same as above but the user doesnt have to provide a file name. It is automatically assumed from Cube's name }
begin
 Result:= AsFasta(ChangeFileExt(FileName, '.FASTA'), Clean);
end;












{--------------------------------------------------------------------------------------------------
                                  PROPERTIES
--------------------------------------------------------------------------------------------------}
function TCubeImport.PropertiesWithName: string;                                                   { inainte sa apelez asta trebuie sa aplic TrimEnds }
begin
 Result:= 'FILE: '+ CRLF
            + '  '+ ShortName+ crlf
            + '  '+ Properties;
end;



function TCubeImport.Properties: string;                                                           { inainte sa apelez asta trebuie sa aplic TrimEnds }
VAR Tags, sComment: string;
    Proc: Double;
begin
 Result:= '';

 if HasChroma
 then
  begin
   Proc:= PercentGC(Bases);

   Result:= Result+ CRLF+ 'GC: '+ Real2Str(Proc, 2)+ '%';
   Result:= Result+ CRLF+ 'AT: '+ Real2Str(100- Proc, 2)+ '%';
   Result:= Result+ CRLF+ ShowBases;

   Result:= Result+ CRLF+ CRLF+ 'QUALITY INFO:';
   if QvComputedInternally                                                                         { If true, it means that the internal algorithm was used to compute the QVs }
   then Result:= Result+ CRLF+ '  Base quality values (QV): computed by '+ AppData.AppName
   else Result:= Result+ CRLF+ '  Base quality values (QV): present.';

   Result:= Result+ CRLF+ '  '+ ShowBasesAfterTrim;

   if QVExist then
    begin
     Result:= Result+ CRLF+ '  Average quality (before trimming): '+ i2s(AverageQVAll);            { Media QV-urilor pt toata Chromatograma (low qv ends included) }
     Result:= Result+ CRLF+ '  '+ ShowEstimatedQuality;                                            { Media QV-urilor pt bucata buna din Chromatograma }
     Result:= Result+ CRLF+ '  Trusted bases (bases with QV over '+ IntToStr(EngTrim1.GoodQVTresh) + '): '+ Real2Str(TrustedBasesPercent)+ '%';
    end;
   Result:= Result+ CRLF+ '  Average peak height: ' + i2s(AveragePeakHeight);
  end
 else
   Result:= Result+ CRLF+ 'No of bases: '+ i2s(NoOfBases);


 { COMMENT }
 if (Comment= '>') OR (Comment= '> ') OR (Comment= '')
 then sComment:= ''
 else sComment:= 'Comment:'+ CRLF+ Comment;

 if CommentOrig <> ''
 then Tags:= 'Abi Tags: '+ CommentOrig;

 { SHOW SHORTER TEXT FIRST }
 if (Tags<> '') OR (sComment<> '') then
   if Length(Tags) > Length(sComment)
   then Result:= Result+ crlf+ sComment+ LBRK+ Tags
   else Result:= Result+ crlf+ Tags   + LBRK+ sComment;

 Result:= System.SysUtils.AdjustLineBreaks(Result, tlbsCRLF)
end;


{ How many trusted bases (bases with QV over GoodQVTresh) are there? }
function TCubeImport.TrustedBasesPercent: Real;
VAR i, TotalGood: Integer;
begin
 TotalGood:= 0;
 for i:= GoodQVStart to GoodQVStart+ NoOfGoodBases-1 DO
  if CellsMX[i].QV >= EngTrim1.GoodQVTresh
  then Inc(TotalGood);
 Result:= ProcentRepresent(TotalGood, NoOfGoodBases);
end;


{ Media QV-urilor pt bucata buna din Chromatograma }
function TCubeImport.ShowEstimatedQuality: string;
begin
 if QVExist
 then
  begin
    Result:= 'Average quality (after trimming): '+ i2s(AverageQV);
    case AverageQV of
     00..14 : Result:= Result+ ' (very poor)';
     15..24 : Result:= Result+ ' (poor)';
     25..39 : Result:= Result+ ' (good)';
     40..59 : Result:= Result+ ' (very good)';
     60..89 : Result:= Result+ ' (excellent)';
    else
     Result:= Result+ ' (incredibly good)';
    end;
  end
 else Result:= 'Sample quality cannot be estimated (base quality info is missing).';
end;



function TCubeImport.ShowQVExist: string;
begin
 if QVExist
 then Result:= 'Quality value info: present'
 else Result:= 'Quality value info: missing!';
end;



function TCubeImport.ShowBases: string;
begin
 Result:= 'Number of bases: '+ i2s(NoOfBases);
end;



function TCubeImport.ShowGoodBases: string;
begin
 Result:= 'Number of good bases: '+ i2s(NoOfBases);
end;


{ Tell how many % of bases are left after TrimEnds }
function TCubeImport.ShowBasesAfterTrim: string;
VAR BasesLeft: Integer;
begin
 if FQVExist then
  begin
   BasesLeft:= round( ProcentRepresent(NoOfGoodBases, NoOfBases) );
   Result:= 'Bases left after end trimming: '+ i2s(NoOfGoodBases)+ ' ('+ i2s(BasesLeft)+ '%)';
  end; 
end;









{--------------------------------------------------------------------------------------------------
   BUILD FROM SCRATCH
   Build a cube from the given parameters.
   After that we have to call ComputeColors.
--------------------------------------------------------------------------------------------------}
procedure TCubeImport.BuildFromScratch;
begin
 Assert(NoOfBases= 0);
 {FILE}
 Comment     := sComments;
 FileName    := sFullName;
 {cube}
 FQVExist    := FALSE;
 NoOfBases   := Length(sBases);
 NoOfSamples := 0;
 setBases(sBases);                                                                                 { This also calls 'DirtyBases:= true' }
 Boot;                                                                                             { This calls DetectVectors }
end;










{--------------------------------------------------------------------------------------------------
   DEBUG
--------------------------------------------------------------------------------------------------}
procedure TCubeImport.DebugToTxt(Path: string);
VAR   Smpl: Integer;
      sFinal, Rw1, Rw2, Rw3, Rw4: string;
begin
 Rw1:= 'Sample:   ';
 Rw2:= 'Bases :   ';
 Rw3:= 'Base pos: ';
 Rw4:= 'CellMX  : '; 

 for Smpl:= ctChromaIndex to NoOfSamples DO
  begin
   Rw1:= Rw1+ i2s(Smpl)+ Tab;
   Rw2:= Rw2+ Char(Sample2Base(Smpl))+ Tab;
   Rw3:= Rw3+ i2s(Chroma[Smpl].Ptr2Base)+ Tab;

   if Chroma[Smpl].Ptr2Base> NoneAssigned
   then Rw4:= Rw4+ Char(Sample2Base(Smpl))+ Tab
   else Rw4:= Rw4+ ' '+ Tab;
  end;

 sFinal:= Rw1+ CRLF+ Rw2+ CRLF+ Rw3+ CRLF+ Rw4
   +CRLF
   +CRLF+ 'Details:'
   +CRLF+ '  NoOfBases '+ i2s(NoOfBases)
   +CRLF+ '  NoOfSamples '+ i2s(NoOfSamples)
   +CRLF+ '  File '+ FileName
   +CRLF;
 StringToFile(Path, sFinal, woOverwrite);
end;



procedure TCubeImport.DebugToBmpSave(SaveTo: string; ShowLinks: Boolean; CONST YDivider: Integer= 10; CONST XMultiplier: Integer= 5);   { creates a bitmap that shows the internal laying of the chromatogram. 10 is a good default value for YDivider. if Path is empty, the routine will autogenerate the path }
VAR bmp: TBitmap;
begin
 if SaveTo= ''
 then SaveTo:= ChangeFileExt(Self.FileName, '.bmp');

 Assert(YDivider> 0);

 bmp:= DebugToBmp(ShowLinks, YDivider, XMultiplier);
 TRY
  bmp.SaveToFile(SaveTo);
 FINALLY
  FreeAndNil(bmp);
 END;
end;



{ Creates a bitmap that shows the internal laying of the chromatogram.
   YDivider is a factor used to lower the height of the samples. 10 is a good default value.
   XMultiplier is a factor to scale the samples on x axis. 4 is a good default value. }

function TCubeImport.DebugToBmp(ShowLinks: Boolean; YDivider, XMultiplier: integer): TBitmap;
CONST BaseLine = 20;                                                                               { linia unde sunt desenate bazele }
VAR   Smpl, SmplX, BaseNr, FGnd, iBase, Highest, iBaseNoGap: Integer;
      BaseDistance: integer;                                                                       { average distance between two bases (in smaples) }
      FirstPointer: Integer;                                                                       { search for the first sample that has a base assigned to it }
      BaseX: Integer;                                                                              { place of the base on the x axis }
begin
 Result:= TBitmap.Create;
 Result.Canvas.Brush.Color:= clBlack;
 Result.Canvas.Brush.Style:= bsClear;
 Result.Canvas.Pen  .Color:= clGray;
 Result.Canvas.Font .Size := 7;
 Result.Width := Self.NoOfSamples * XMultiplier+ 20;                                               { mai las 20 de pixeli la coada casa vad si ultima baza, altfel ea e desenata afara din ecran }
 Result.Height:= 300;

 FGnd := Result.Height- 14;
 BaseDistance:= 12* XMultiplier;                                                                   { average distance between two bases (in smaples) }

 FirstPointer:= 0;
 for Smpl:= ctChromaIndex to NoOfSamples DO
   if Chroma[Smpl].Ptr2Base> NoneAssigned then
    begin
     FirstPointer:= Smpl* XMultiplier;                                                             { search for the first sample that has a base assigned to it }
     Break;
    end;

 { LABELS }
 Result.Canvas.Font.Color:= clWhite;
 Result.Canvas.TextOut(3, 75, 'Y scale divided by '+ i2s(YDivider));
 Result.Canvas.TextOut(3, 87, 'X scale multiply by '+ i2s(XMultiplier));
 Result.Canvas.Font.Color:= clSilver;
 Result.Canvas.TextOut(3, 99, 'S = sample no. ');
 Result.Canvas.TextOut(3, 111, 'P = ptr to base. ');

 Result.Canvas.Font.Color:= clOrange;
 Result.Canvas.TextOut(0, BaseLine+25, 'ChrmMX:');                                                 { Label for bases from CHROMAMX }
 Result.Canvas.Font.Color:= clLime;
 Result.Canvas.TextOut(0, BaseLine+5, 'CellMX:');                                                  { Label for bases from CELLMX }

 with Result.Canvas DO
 if ShowLinks then
  for Smpl:= ctChromaIndex to NoOfSamples DO
   begin
    Highest:= 0;
    SmplX:= Smpl* XMultiplier;                                                                      { place of the sample on the x axis }

    { highest peaks from all 4 (ACTG) }
    if Chroma[Smpl].HeightA> Highest then Highest:= Chroma[Smpl].HeightA;
    if Chroma[Smpl].HeightC> Highest then Highest:= Chroma[Smpl].HeightC;
    if Chroma[Smpl].HeightG> Highest then Highest:= Chroma[Smpl].HeightG;
    if Chroma[Smpl].HeightT> Highest then Highest:= Chroma[Smpl].HeightT;
    Highest:= FGnd -(Highest div YDivider);

    { Show links between highest peak and its associated base }
    if Chroma[Smpl].Ptr2Base> NoneAssigned then
     begin
      BaseX:= FirstPointer+ Chroma[Smpl].Ptr2Base* BaseDistance;                                    { BaseX -> place of the base on the x axis }
      Pen.Color:= rgb(45, 45, 55);
      MoveTo(SmplX, Highest-2);                                                                    { from peak }
      LineTo(BaseX, BaseLine+ 54);                                                                 { to base  (+54 e ca sa fac sageata sa se termine un pic sub baze) }
    end;
   end;

 BaseNr:= 0;
 with Result.Canvas DO
 for Smpl:= ctChromaIndex to NoOfSamples DO
  begin
   Highest:= 0;
   SmplX:= Smpl* XMultiplier;                                                                      { place of the sample on the x axis }

   { draw samples }
   Pixels[SmplX, FGnd- Chroma[Smpl].HeightA div YDivider]:= clGreenDark;
   Pixels[SmplX, FGnd- Chroma[Smpl].HeightC div YDivider]:= clRed;
   Pixels[SmplX, FGnd- Chroma[Smpl].HeightG div YDivider]:= clBlue;
   Pixels[SmplX, FGnd- Chroma[Smpl].HeightT div YDivider]:= clSilver;

   { highest peaks from all 4 (ACTG) }
   if Chroma[Smpl].HeightA> Highest then Highest:= Chroma[Smpl].HeightA;
   if Chroma[Smpl].HeightC> Highest then Highest:= Chroma[Smpl].HeightC;
   if Chroma[Smpl].HeightG> Highest then Highest:= Chroma[Smpl].HeightG;
   if Chroma[Smpl].HeightT> Highest then Highest:= Chroma[Smpl].HeightT;
   Highest:= FGnd -(Highest div YDivider);

   { draw pointers from peaks to 'bases line' }
   if Chroma[Smpl].Ptr2Base> NoneAssigned then
    begin
     inc(BaseNr);
     BaseX:= FirstPointer+ Chroma[Smpl].Ptr2Base* BaseDistance;                                    { BaseX -> place of the base on the x axis }

     Result.Canvas.Font.Color:= clWhite;

     { Peak's associated base }
     TextOut(SmplX, Highest-17, Char(Sample2Base(Smpl)));
     { Peak's position }
     TextOut(SmplX, Highest-32, 'S:'+i2s(Smpl));
     { Peak's pointer }
     TextOut(SmplX, Highest-47, 'P:'+ i2s(Chroma[Smpl].Ptr2Base));

     { BASES - Draw base from CHROMAMX }
     Result.Canvas.Font.Color:= clOrange;
     TextOut(BaseX-3, BaseLine+25, Char(Sample2Base(Smpl)));
     { Draw base's number }
     Result.Canvas.Font.Color:= clPurpleLight;
     TextOut(BaseX-3, BaseLine+38, i2s(BaseNr));
   end;
  end;

 {= Draw base from CELLMX =}
 iBase:= 0;
 iBaseNoGap:= 0;
 BaseX:= 1;
 with Result.Canvas DO
  while iBase< NoOfBases DO
   begin
     inc(iBase);
     if Base[iBase]<> Gap
     then
      begin
       inc(iBaseNoGap);
       BaseX:= FirstPointer+ iBaseNoGap* BaseDistance;                                             { BaseX -> place of the base on the x axis }
       Result.Canvas.Font.Color:= clGreen;
       TextOut(BaseX-3, BaseLine+5, Char(Base[iBase]));
       Result.Canvas.Font.Color:= clLime;
       TextOut(BaseX-3, BaseLine-8, i2s(iBase));
      end
     else
      begin
       BaseX:= BaseX+  12*XMultiplier;                                                             { inghesuie Gap-urile }
       Result.Canvas.Font.Color:= clRed;
       TextOut(BaseX-3, BaseLine+5, Char(Base[iBase]));
       Result.Canvas.Font.Color:= clRed;
       TextOut(BaseX-3, BaseLine-8, i2s(iBase));
      end;

   end;
end;



{SAVE cube BITMAP}
procedure TCubeImport.SaveChromaAsBmp(const nume: string);
Var i: Integer;
    img: TBitmap;
CONST
   GND= 230;
begin
 img:= TBitmap.Create;
 Img.Width:= NoOfSamples;
 Img.Height:= 300;

 with Img.Canvas DO
  begin
   Font.Size:= 6;
   Font.Name:= 'Arial';
   MoveTo(0, GND); Pen.Color:= clGreen;
   for i:= ctChromaIndex to NoOfSamples DO
     LineTo(i, GND- Chroma[i].HeightA DIV 30);
   MoveTo(0, GND); Pen.Color:= clBlue;
   for i:= ctChromaIndex to NoOfSamples DO
     LineTo(i, GND- Chroma[i].HeightC DIV 30);
   MoveTo(0, GND); Pen.Color:= clBlack;
   for i:= ctChromaIndex to NoOfSamples DO
     LineTo(i, GND- Chroma[i].HeightG DIV 30);
   MoveTo(0, GND); Pen.Color:= clRed;
   for i:= ctChromaIndex to NoOfSamples DO
     LineTo(i, GND- Chroma[i].HeightT DIV 30);

   for i:= ctChromaIndex to NoOfSamples DO
    if Chroma[i].Ptr2Base> NoneAssigned then
      begin
       TextOut(i, GND+20, Char(Sample2Base(i)));
       TextOut(i, GND+35, i2s(Chroma[i].Ptr2Base));
      end;
  end;

 img.SaveToFile(nume);
 FreeAndNil(img);
end;



procedure TCubeImport.FakeChroma;                                                                  { FAKE A CHROMA  -  PeakHeight should be about set at 700. ATENTIE: trebuie sa apelez "Randomize" inainte a a apela "FakeChroma" }
VAR Sample, CurBase: Integer;                                                                      { Issues if I open multiple files. }
    Baza: TBase;
    Sinoid1, Sinoid2, Sinoid3, Sinoid4, Sinoid5: Integer;
    sn1, sn2, sn3, sn4, sn5: Integer;
    Ran: integer;    
const
    BaseDistance= 12;
begin
 Clear;
 Sinoid1:= round(PeakHeight* sin(1/5));                                                            { Tabela de constante }
 Sinoid2:= round(PeakHeight* sin(2/5));
 Sinoid3:= round(PeakHeight* sin(3/5));
 Sinoid4:= round(PeakHeight* sin(4/5));
 Sinoid5:= round(PeakHeight* sin(1));                                                              {5/5}

 sn1:= Sinoid1;
 sn2:= Sinoid2;
 sn3:= Sinoid3;
 sn4:= Sinoid4;
 sn5:= Sinoid5;

 { Mark all bases as non-gray }
 for CurBase:= 1 TO NoOfBases DO
  CellTrimmed[CurBase]:= FALSE;

 if LastGoodBase<= GoodQVStart
 then LastGoodBase:= NoOfBases;

 {BUILD CHROMA}
 Sample:= 0;
 NoOfSamples:= (NoOfBases)* BaseDistance;                                                          { 1 Rezerv memorie pentru ele. 2 Aflu numarul teoretic de sample-uri din care scad GAP-urile caci pentru GAPuri nu am PEAK-uri }
 FQVExist:= TRUE;
 Assert(length(BASES)= NoOfBases, 'Length(BASES)=NoOfBases        FAILED');

 for CurBase:= 1 to NoOfBases DO
  begin
   inc(Sample, BaseDistance);

   { Fake the QV }
   if  (CurBase> GoodQVStart)
   AND (CurBase< LastGoodBase)
   then CellQV[CurBase]:= 30                                                                       { BIG QV FOR MIDLE }
   else CellQV[CurBase]:= 10;                                                                      { LOW QV for heads}

   Chroma[Sample].Ptr2Base:= CurBase;

   if Randomiz then
    begin
     Ran:= Random(200);
     if NOT ( (CurBase> GoodQVStart) AND (CurBase< LastGoodBase) )
     then Ran:= -Ran;

     sn1:= Sinoid1+ Ran;
     sn2:= Sinoid2+ Ran;
     sn3:= Sinoid3+ Ran;
     sn4:= Sinoid4+ Ran;
     sn5:= Sinoid5+ Ran;
    end;

   Baza:= BASES[CurBase];   
   if Baza= 'A' then                                                                               { pntru toate bazele ACGT desenez un frumos sinus }
     begin
      Chroma[Sample-4].HeightA:= sn1;
      Chroma[Sample-3].HeightA:= sn2;
      Chroma[Sample-2].HeightA:= sn3;
      Chroma[Sample-1].HeightA:= sn4;
      Chroma[ Sample ].HeightA:= sn5;
      Chroma[Sample+1].HeightA:= sn5;
      Chroma[Sample+2].HeightA:= sn4;
      Chroma[Sample+3].HeightA:= sn3;
      Chroma[Sample+4].HeightA:= sn2;
      Chroma[Sample+5].HeightA:= sn1;
     end else
   if Baza= 'C' then
     begin
      Chroma[Sample-4].HeightC:= sn1;
      Chroma[Sample-3].HeightC:= sn2;
      Chroma[Sample-2].HeightC:= sn3;
      Chroma[Sample-1].HeightC:= sn4;
      Chroma[ Sample ].HeightC:= sn5;
      Chroma[Sample+1].HeightC:= sn5;
      Chroma[Sample+2].HeightC:= sn4;
      Chroma[Sample+3].HeightC:= sn3;
      Chroma[Sample+4].HeightC:= sn2;
      Chroma[Sample+5].HeightC:= sn1;
     end else
   if Baza= 'T' then
     begin
      Chroma[Sample-4].HeightT:= sn1;
      Chroma[Sample-3].HeightT:= sn2;
      Chroma[Sample-2].HeightT:= sn3;
      Chroma[Sample-1].HeightT:= sn4;
      Chroma[ Sample ].HeightT:= sn5;
      Chroma[Sample+1].HeightT:= sn5;
      Chroma[Sample+2].HeightT:= sn4;
      Chroma[Sample+3].HeightT:= sn3;
      Chroma[Sample+4].HeightT:= sn2;
      Chroma[Sample+5].HeightT:= sn1;
     end else
   if (Baza= 'N') OR (Baza= Gap) then     {OR (Baza= Gap Asm)}                          { pentru GAP-uri pun peak-uri foarte proaste }
     begin
      Chroma[Sample-4].HeightT:= 10;
      Chroma[Sample-3].HeightT:= 20;
      Chroma[Sample-2].HeightT:= 30;
      Chroma[Sample-1].HeightT:= 40;
      Chroma[ Sample ].HeightT:= 50;
      Chroma[Sample+1].HeightT:= 50;
      Chroma[Sample+2].HeightT:= 40;
      Chroma[Sample+3].HeightT:= 30;
      Chroma[Sample+4].HeightT:= 20;
      Chroma[Sample+5].HeightT:= 10;
     end else
   { ORICE ALTA LITERA, INTRA LA G - AR TREBUIE SA FOLOSESC AMBIGUITY CODE! }
     begin
      Chroma[Sample-4].HeightG:= sn1;
      Chroma[Sample-3].HeightG:= sn2;
      Chroma[Sample-2].HeightG:= sn3;
      Chroma[Sample-1].HeightG:= sn4;
      Chroma[ Sample ].HeightG:= sn5;
      Chroma[Sample+1].HeightG:= sn5;
      Chroma[Sample+2].HeightG:= sn4;
      Chroma[Sample+3].HeightG:= sn3;
      Chroma[Sample+4].HeightG:= sn2;
      Chroma[Sample+5].HeightG:= sn1;
     end;
  END;                                                                                             {WHILE}
end;




END.
