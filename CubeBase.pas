
UNIT CubeBase;

{==================================================================================================
 Heracle BioSoft SRL
 2016.11.01

 Object capabilities:
    + Soft editing capabilities (Replace base, Add, Delete base )
    + User editing capabilities (Replace base)
    + Search
    + Vectors (GoodBasesNoVector, Detect vectors, List of detected vectors )

 Documentation:
    * SizeOf the RSample structure is 16 bytes because the structure is
      aligned to 8 bytes. If I use packed then the structure will be 12 bytes.

==================================================================================================}

INTERFACE

USES
   System.SysUtils, System.StrUtils, Vcl.Graphics, ccCore, ccINIFile, CubicDNA, ReadBasicSample, clRamLog, ccRichLog;

CONST
   NoneAssigned  = 0;                                                                              { Unde puncteaza pointerul Chroma.Ptr2Base cand nici o baza nu ii e asignata }
   ctChromaIndex = 1;
   ChromaIndex= ctChromaIndex;                                                                     { Chroma MX e indexata in 1 (dar pe viitor am de gand sa o indexez in zero) }
   ctCellsIndex  = 1; CellsIndex= ctCellsIndex;                                                    { Cells  MX e indexata in 1 }

TYPE
 TChanged    = procedure (Sender: TObject) of object;
 TMismStat   = (msNormal,                                                                          { change its name to msNone                                            { No mismatch, mutation or SNP detected at this cell }
                msMisUnsolved,                                                                     { Mismatch detected. The user took no action to solve it }
                msMisSolved);                                                                      { Mismatch detected. The user approved the base suggested by DNA Baser or entered his own choice }                                                                             { This is valid for contig only!!!! One of the sequence has a base that the other is missing. }
 TCtgMutation= (cmNone,                                                                            { No mutation in this cell }
                cmDifference,                                                                      { This is valid for contig only!!!! A mutation was detected here (two sequences in assembly grid have different bases at this point) }
                cmIndel);                                                                          { This is valid for contig only!!!! One of the sequence has a base that the other is missing. }

 {IMPORTANT DON'T CHANGE THE SIZE AND ORDER OF THESE RECORDS }

 RBaseProb = record
    A: BYTE;                                                                                       { Probability of it being an A }
    C: BYTE;                                                                                       { Probability of it being an C }
    G: BYTE;                                                                                       { Probability of it being an G }
    T: BYTE;                                                                                       { Probability of it being an T }
  end;

{
There are 3 cases in which a base can differe from its original (ABI/SCF) value.
 1. It was recalculated by DNA Baser base caller (the N bases)
 2. It was recalculated by DNA Baser base caller (any base that was now recomputed to a SNP )
 3. It was edited by the user in Grid
 }

 RCell5 = record                                                                                   { Current version (v4.1) }
  private
    procedure Reverse;                                                                             { Reverse the base and its probability. To be used only in Cube.Reverse }
  public
    Base        : TBase;                                                                           { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
    QV          : Byte;
    MismStatus  : TMismStat;                                                                       { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
    Trimmed     : Boolean;

    EdtOrig     : TBase;                                                                           { Base before user edit. Makes sense for input samples only, not also for contig }
    EditPtr     : Integer;                                                                         { Makes sense for AsmJob. Points to the EditMX entry coresponding to the base edit (if this base was indeed edited) }
    SnpOrig_    : TBase;   { UNUSED }                                                             { This field stores the original base as imported from SCF/ABI file. If different than NoBase then it means that a SNP was detected here.   DEL: TConsensor.Conses will set this field to the same value as 'Base'  }
    CellProb    : RBaseProb;                                                                       { Probability of it being an A/C/G/T }

    Bookmark    : Boolean;
    Vector      : Boolean;                                                                         { True if this cell has a Vector assigned to it }
    ColorB      : TColor;                                                                          { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    ColorT      : TColor;                                                                          { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    Ptr2Smpl    : Integer;                                                                         { position of the sample associated with this base }
    CtgMutation : TCtgMutation;
    SnpAreaRatio: Byte;                                                                            { Mutant peak to reference peak area ratio (%). I need it for the Log }
    SnpOverlap  : Byte;                                                                            { Overlap between mutant peak and reference peak (%). I need it for the Log  }
    procedure Clear;
    function BaseOrig: TBase;                                                                      { If different than NoBase then it means that a SNP was detected here. It does not store the SNP base but a backup of the original base as imported from SCF/ABI file }
  end;

 RSample = record
    HeightA    : Word;                                                                             { Amplitudinea esantionelor }
    HeightC    : Word;                                                                             { Atentie! Sample ACGT erau Smallint insa l-am trecut la Word din cauza ca un sample nu are cum sa fie negativ. }
    HeightG    : Word;
    HeightT    : Word;
    Ptr2Base   : Integer;                                                                          { Pointer catre o baza. S-ar putea sa fie nul sau s-ar putea chiar sa puncteze catre o baza. }
  end;

 ACells     = ARRAY of RCell5;                                                                     { for DNA_Baser v4.1 cubes (current version) - Baser v4.9 and up  }
 AChromaMX  = ARRAY of RSample;                                                                    { toate cele 4 trace-uri suprapuse formeaza SAMPLE MATRIX (care e ca un BITMAP) }
 XTrace     = ARRAY of Word;

 RTrace     = record
   Height: XTrace;
   Color : TBase;
 end;

 TCubeAbstract = class(TBasicSample)
  private
   FSamples  : Integer;
   FNoOfBases: Integer;                                                                            { CellMX e indexat in 1 nu in 0 }       { Number of bases in Bases matrix.  IMPORTANT. BasesNr e diferit de BasesNrOrig pentru ca userul/asamblarea poate adauga baze noi in plus. Cand import un fisier o sa folosesc BasesNrOrig, insa dupa aceea o sa fol. BasesNr }
   FChanged  : TChanged;
   FReversed : Boolean;                                                                            { True arata ca i s-a aplicat functia reverse }
   function  getBases          : BaseString;
   function  getBase           (Index: integer): TBase;                                            { Base }
   procedure setBase           (Index: integer; Value: TBase);                                     { Set the base, call Changed which rebuilds the 'Bases' string }
   function  getEdtOrig        (Index: integer): TBase;
   procedure setEdtOrig        (Index: integer; Value: TBase);
   function  getEditPoint      (Index: integer): Integer;
   procedure setEditPoint      (Index, Value: Integer);
   function  getSample         (Index: integer): RSample;
   procedure setSample         (Index: Integer; CONST Value: RSample);
   function  getPtr2Smpl       (Index: integer): Integer;
   procedure setPtr2Smpl       (Index: integer; Value: Integer);
   function  getHasVect        (Index: integer): Boolean;
   procedure setHasVect        (Index: integer; Value: Boolean);
   function  getColorB         (Index: integer): Tcolor;
   procedure setColorB         (Index: integer; Value: Tcolor);
   function  getColorT         (Index: integer): Tcolor;
   procedure setColorT         (Index: integer; Value: Tcolor);
   function  getBookmark       (Index: integer): Boolean;
   procedure setBookmark       (Index: integer; Value: Boolean);
   procedure setNrSamples      (CONST Value: Integer);

   { GOOD QV }
   function  getQV             (Index: integer): byte;
   procedure setQV             (Index: integer; Value: byte);                                      { QV }
   function  getCellProb       (Index: integer): RBaseProb;
   function  getBaseGr         (Index: integer): Boolean;
   procedure setBaseGr         (Index: Integer; Value: Boolean);
   function  getMismatch       (Index: integer): TMismStat;
   procedure setMismatch       (Index: integer; Value: TMismStat); virtual;
   procedure setLastGoodBase   (LastGoodBase: Integer);
   procedure setGoodQVStart    (Value: Integer);
   function  getGoodBases      : BaseString;
   function  getGoodQVEnd      : Integer;
   function  getGoodQVStart    : Integer;
   function  getGoodStSmpl     : Integer;                                                                  { same as GoodQVStart but returns a sample nr instead of base nr }
   function  getGoodEnSmpl     : Integer;
   {}
   procedure buildGoodBases;
  protected
   CellsMX          : ACells;                                                                              { CellsMX E INDEXATA IN 1 }
   FGoodBases       : BaseString;                                                                          { <- valoarea este refacuta ori de cate ori bazale au fost alterate/ modificate }
   DirtyBases       : Boolean;
   DirtyGoodBases   : Boolean;
   procedure setNoBases        (CONST Value: Integer);        override;
   function  getNoBases: Integer;                             override;
   procedure setBases          (CONST Bases: BaseString);                                                  {NOT USED} { importa un sir de baze. util cand import un obiect de tip SEQ/FASTA care nu are chromatograma si la care GenereazaCellMx nu merge  }
   procedure SetReciprocPointers(BaseIndex, SampleIndex: Integer);
   procedure CellMxReset;
   function  SampleHasBaseAssignedArround(CONST Sample: Integer): Boolean;                                 { I check few positions arround this sample because for some strange reasons, in ABI, the peaks are shifted to the left with 1 or 2 pixels }
   function  DistanceBetween   (CONST Base1, Base2: Integer): Integer;                                     { The distance (in samples) between two bases }
   function  GetMiddlePoint    (CONST FirstBase: Integer; OUT MustCreate: Boolean): Integer;
  public
   Chroma          : AChromaMX;                                                                            { E INDEXATA IN 1 }
   VectorColor     : TColor;
   VectorShow      : Boolean;                                                                              { Show Vectors (in blue color) or not }
   constructor Create(aLog: TRamLog);
   procedure Clear; override;
   procedure buildBases;
   {VECTORS}
   procedure DetectVectors;     override;
   procedure VectorsClearColor; override;
   function  GoodBasesNoVector: BaseString;
   function  VectorAtBase(BasePos: Integer): string;                                                       { Returns the name of the vector associated with the specified base }
   function  LeftVectorExists: Boolean;
   function  RightVectorExists: Boolean;
   function  BaseIsStriked(BasePos: Integer): Boolean;                                                     { Returns true if the specified base will be cut from contig because it is located in vector or previous to the vector }
   procedure CleanFromTo (OUT FromBase, ToBase: Integer);                                                  { Shows from which to which base I have to cut in order to remove the vectors and the low quality ends }

   {CELLS MX}
   function  Cell              (CONST Index: integer): RCell5;
   procedure CellSet           (CONST Index: Integer; CONST NewValue: RCell5);
   function  CellProperties    (CONST BasePos: Integer): string;
   procedure CellReset         (CONST BasePos: integer);
   property  CellQV            [Index: Integer]: Byte      read getQV        write setQV;
   property  CellProbab        [Index: Integer]: RBaseProb read getCellProb;
   property  CellTrimmed       [Index: Integer]: Boolean   read getBaseGr    write setBaseGr;              { Returns true if the base is gray }
   property  CellEdtOrig       [Index: Integer]: TBase     read getEdtOrig   write setEdtOrig;
   property  CellEditPtr       [Index: Integer]: Integer   read getEditPoint write setEditPoint;           { old name: CellPtr2EditMX}  { Points to the EditMX entry coresponding to the base edit (if this base was indeed edited) }
   property  CellMism          [Index: Integer]: TMismStat read getMismatch  write setMismatch;            { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
   property  CellColorB        [Index: Integer]: Tcolor    read getColorB    write setColorB;              { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
   property  CellColorT        [Index: Integer]: Tcolor    read getColorT    write setColorT;              { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
   property  CellHasVector     [Index: Integer]: Boolean   read getHasVect   write setHasVect;             { True if this cell has a Vector assigned to it }
   property  CellBookmark      [Index: Integer]: Boolean   read getBookmark  write setBookmark;

   {BASES}
   property  Base              [Index: Integer]: TBase     read getBase      write setBase;                { Indexed in 1 }
   property  Base2Smpl         [BasePos: Integer]: Integer read getPtr2Smpl  write setPtr2Smpl;            { Pointer to sample - Position of the sample associated with this base }   { intoarce sample-ul corespunzator pt baza data ca parametru }
   property  Bases             : BaseString                read getBases;                                  { <- valoarea este refacuta ori de cate ori bazale au fost alterate/ modificate }

   { MUTANT BASE }
   function  CellSnpOverlap    (Index: Integer): Byte;                                                     { Overlap between mutant peak and reference peak }
   function  CellSnpArea       (Index: Integer): Byte;                                                     { Mutant peak to reference peak area ratio. I need it for the Log }
   function  CellIsIUPAC        (Index: integer): Boolean;                                                { Shows if the base at the specified pos is an IUPAC base }
   function  CellIsIUPACDetected(Index: Integer): Boolean;                                                { Shows if the base at the specified pos was detected BY DNA BASER as a SNP base }

   function  CellIsMutant      (Index: integer): Boolean;                                                  { Returns TRUE if the cell has ANY kind of mutation (SNP, Indel, Difference }
   function  CellIsMutDiff     (Index: integer): Boolean;                                                  { Valid only for a contig base! Returns true if the contig base is an IUPAC base bacause the input sequences are different at this position }
   function  CellIsMutIndel    (Index: integer): Boolean;

   {BASES-EDIT}
   procedure ReinsertGoodBases (CONST NewMiddle: BaseString);
   procedure ReplaceBase       (CONST BasePos: integer; CONST NewBase: TBase; CONST QV: Byte);             { replace a base in chroma, in UnitsMX si BASES }
   procedure ReplaceBase2Gap   (CONST BasePos: integer);                                                   { replace the specified base with a gap }
   procedure SwapCells         (CONST Cell1, Cell2: Integer);
   procedure AppendBaseAtBeg   (CONST NewBase: TBase; CONST QV: Byte);                                     { Append base at the beginning of the chromatogram }
   procedure AppendBaseAtEnd2  (CONST NewBase: TBase; CONST QV: Byte); overload;                           { Append base at the end of the chromatogram }
   procedure AppendCellAtEnd   (CONST NewCell: RCell5);                 overload;                          { Append base at the end of the chromatogram }
   procedure InsertBase        (CONST NewBase: TBase; CONST QV: byte; CONST BeforeBase: Integer);          { The base at position 'BeforeBase' and all bases after this position will be pushed to the right. The empty placed created at 'BeforeBase' position will be filled by the new base }

   procedure TruncateTo        (CONST BasePos: Integer);                                                   { Deletes are bases/samples from the begining of the sample to BasePos (including). The chromatogram will be shortened}
   procedure TruncateFrom      (CONST BasePos: Integer);
   procedure RemoveBase        (CONST BasePos: Integer);                                                   { BaseNr - arata baza care va fi stearsa }
   function  RemoveBasePossible(CONST BasePos: Integer): boolean;                                          { is it possible to delete this base? }
   procedure RemoveAllGaps;

   {BASES-MANUAL EDITS}
   procedure EditBase          (CONST BasePos: integer; NewBase: TBase; CONST EditMxIndex: Integer);                       overload;      { Set the base and then set the Edited to TRUE. Called when the user manually edits a base in TContigGrid }
   procedure EditBase          (CONST BasePos: integer; NewBase: TBase; CONST EditMxIndex: Integer; Mismatch: TMismStat);  overload;

   {BASES-GOOD BASES}
   property  GoodBases         : BaseString   read getGoodBases;                                           { <- valoarea este refacuta ori de cate ori bazele au fost alterate/ modificate sau s-a aplicat rimming engine }
   function  NoOfGoodBases     : integer;
   function  GoodBasesRC       : BaseString;                                                               { RC of GoodBases string }
   property  GoodQVStart       : integer  read getGoodQVStart  write setGoodQVStart;                       { including     |??????|kkkkkkkkkkkkkkkk|???????|      } { ATENTIE!  E indexat in 1 }  { The place where the first good base is located }
   property  LastGoodBase      : integer  read getGoodQVEnd    write setLastGoodBase;                      { including            |G1            G2|              } { Original name was GoodQVEnd }
   property  GoodQVStartS      : integer  read getGoodStSmpl;                                              { same as GoodQVStart but returns a sample nr instead of base nr }
   property  GoodQVEndS        : integer  read getGoodEnSmpl;                                              { same as GoodQVEnd but returns a sample nr instead of base nr }
   function  HasAmbiguousBases : Boolean;                                                                  { True if there are non ACGT bases }

   {QUALITY}
   function  AveragePeakHeight : Integer;                                                                  { Media peak-urilor pt toata Chromatograma }
   function  GoodBasesPercent  : Real;                                                                     { Cat % din baze sunt bune (GoodQV) dupa TrimEngine }

   {CONVERSION}
   function  Sample2QV         (CONST Sample: integer): byte;                                              { intoarce QV-ul corespunzator pt sample-ul dat ca parametru }
   function  Sample2Base       (CONST Sample: integer): TBase;                                              { intoarce baza corespunzatoare pt sample-ul data ca parametru }  { Old name: BaseAtSample}
   function  Sample2BaseI      (CONST Sample: integer): Integer;
   function  Sample2BaseApprox (CONST Sample: integer): integer;                                           { Base pos corresponding to this sample pos. If there is no base assigned to this sample, than find the closes neigbor }

   {SAMPLES/CHROMA}
   function  SampleHasBaseAssigned  (CONST Sample: integer): Boolean;                         overload;     { This checks if the specified sample has a base assigned to it (doesn't matter which color is that base) }
   function  SampleHasBaseAssigned  (CONST Sample: integer; CONST BaseColor: TBase): Boolean; overload;     { This checks if the specified sample has a base (of the specified color) assigned to it }
   function  SampleHasBaseAssignedEx(CONST Sample: Integer; CONST BaseColor: TBase): Boolean;               { Same as above but it ignores the BaseColor parameter to Ns. In other words it always returns TRUE if there is a base at this sample and that base is N. }
   function  HasChroma         : Boolean;
   property  NoOfSamples       : Integer                  read FSamples  write setNrSamples;               { Number of elements in Samples matrix. Indexata in 1 }
   property  Sample            [Index: Integer]: RSample  read getSample write setSample;                  { Access to chromatogram's samples }
   function  Base2SampleHeight (CONST BaseNo: Integer)         : Integer;   overload;                      { Gets a base as input and returns the sample height for the peak associated to that base }
   function  Base2SampleHeight (CONST BaseNo, aSample: Integer): Integer;   overload;
   function  Base2SampleHeight (CONST aBase: TBase; aSample: Integer): Integer; overload;                  { Gets a base as input and returns the sample height for the peak associated to that base }
   function  HighestTrace      (aSample: Integer): Integer;                                                 { Check the ACGT traces for the highest value }
   procedure ExtractTrace      (TraceColor: TBase; OUT Trace: RTrace);
   {}
   property  Reversed          : Boolean read FReversed write FReversed;                                   { True if the sequence was complement-reversed }
   procedure Reverse;                                                                                      { Asta apeleaza la urma pe Changed si ApplyTrimEng }
   procedure Ulimination;                                                                                  { Replaces the 'U' bases with 'T' }
   function  Search            (CONST SearchStr: BaseString; Offset: Integer): Integer; override;
   property  OnChange          : TChanged  read FChanged  write FChanged;
 end;



IMPLEMENTATION
USES
   cmMath, VectorProcessing;






{===============================================================================
   CREATE
===============================================================================}
constructor TCubeAbstract.Create(aLog: TRamLog);                                                                  { Unless you are careful, the object might not be fully constructed when the method is called. To avoid any problems, you should override the AfterConstruction method and use that for any code that needs to wait until the object is fully constructed. If you override AfterConstruction, be sure to call the inherited method, Top. }
begin
 inherited Create(aLog);

 DirtyBases      := TRUE;                                                                          { DirtyBases= True, implicitelly means DirtyGoodBases= true }
 DirtyGoodBases  := TRUE;
 FGoodBases      := '';
 FNoOfBases      := 0;
 NoOfSamples     := 0;
 FReversed       := FALSE;

 { VECTORS }
 VectorShow  := TRUE;                                                                              { show Vectors in blue }
 VectorColor := CellColorVector;
end;


{ 'Clear' is called when:
     * creez un cub abstract (TCubeBase, TCubeImport, TCubObj), atunci apelez Clear imediat dupa Create.
     * creez un cubpersistent (TCubObjEx), atunci apelez Clear cand dau LoadFromFile sau Assign (daca o sa implementez vreodata Assign). }
procedure TCubeAbstract.Clear;
begin
 inherited;

 DirtyBases      := TRUE;                                                                          { DirtyBases= True, implicitelly means DirtyGoodBases= true }
 DirtyGoodBases  := TRUE;
 FGoodBases      := '';
 FNoOfBases      := 0;
 SetLength(CellsMX, 0);
 NoOfSamples     := 0;
 FReversed       := FALSE;
end;








{--------------------------------------------------------------------------------------------------
   NR OF BASES
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.getNoBases: Integer;
begin
 Result:= FNoOfBases;
end;

procedure TCubeAbstract.setNoBases(CONST Value: Integer);
begin
 Assert(Value> 0, '(setNoBases) Invalid number of bases: '+ i2s(Value));

 if FNoOfBases<> value then
  begin
   FNoOfBases:= Value;
   SetLength(CellsMX, FNoOfBases+ ctCellsIndex);   { For a dynamic arrays, SetLength reallocates the array referenced by S to the given length. Existing elements in the array are preserved and newly allocated space is set to 0 or NIL. }
  end;                                                                                             { procedura asta are grija sa adune un 1 pt ca CellMX e indexat in 1 nu in 0 }
end;


procedure TCubeAbstract.setBases(CONST Bases: BaseString);
VAR i: Integer;
begin
 { NOT USED }
 Assert(Length(Bases)= NoOfBases, 'TCubeAbstract.SetBases');
 Assert(Length(Bases)> 0, 'TCubeAbstract.SetBases');
 for i:= 1 to Length(Bases)
  DO CellsMX[i].Base:= Bases[i];

 DirtyBases:= TRUE;
end;




{--------------------------------------------------------------------------------------------------
   DIRTY - BASES ARE DIRTY
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.getBases: BaseString;
begin
 if DirtyBases
 then buildBases;

 Result:= FBases;
end;


function TCubeAbstract.getGoodBases: BaseString;
begin
 if DirtyBases
 then buildBases;                                                                                  { Why is this necessary: The FASTA files (non-qv) the GoodBases is equal with Bases. However, in some places I might call directly GoodBases without building Bases fists. }

 if DirtyGoodBases
 then buildGoodBases;                                                                               { daca nu am apelat inca ApplyTrimEngine, GoodBases s-ar putea sa fie gol }

 Result:= FGoodBases;
end;




{--------------------------------------------------------------------------------------------------
   DIRTY - BUILD BASES
--------------------------------------------------------------------------------------------------}
procedure TCubeAbstract.buildBases;
VAR i: Integer;
begin
 Assert(NoOfBases> 0, 'TCubeAbstract.NoOfBases < 0 ');

 FBases:= '';
 for i:= ctCellsIndex to NoOfBases
  DO FBases:= FBases+ CellsMX[i].Base;

 DirtyBases:= FALSE;
 buildGoodBases;

 if Assigned(FChanged)
 then FChanged(Self);                                                                              { Invalidez randul. da refresh ca sa apara Vectorul cu albastru }
end;


procedure TCubeAbstract.buildGoodBases;                                                            { In implementarea curenta, un fisier poate sa aiba GoodQVStart> 0 chiar daca chromatograma nu contine QV-uri. Asta e posobil din cauza ca userul a marcat MANUAL bazele cu low quality }
VAR i: Integer;
begin
 FGoodBases:= '';
 for i:= GoodQVStart to LastGoodBase DO                                                            { ...calculeaza GoodBases }
  begin
   FGoodBases:= FGoodBases+ Base[i];
   Assert(NOT CellTrimmed[i], '[Error in TCubeAbstract.buildGoodBases]'+CRLF+'Bad QV found instead of a good QV. CRLF Place: CubObj.GoodBases');
  end;

 DirtyGoodBases:= FALSE;
end;




{--------------------------------------------------------------------------------------------------
   BASE ACCESS
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.getBase(Index: integer): TBase;
begin
 Assert(Index>= ctCellsIndex, 'TCubeAbstract.getBase - Index:'+ i2s(Index)+ '  NoOfBas:'+ i2s(NoOfBases));
 Assert(Index<= NoOfBases,    'TCubeAbstract.getBase - Index:'+ i2s(Index)+ '  NoOfBas:'+ i2s(NoOfBases));
 Result:= CellsMX[Index].Base;
end;


procedure TCubeAbstract.setBase(Index: integer; Value: TBase);                                     { Set the base then call Changed which rebuilds the 'Bases' string }
begin
 Assert(Index>= ctCellsIndex, 'TCubeAbstract.setBase - Index:'+ i2s(Index)+ '  NoOfBas:'+ i2s(NoOfBases));
 Assert(Index<= NoOfBases,    'TCubeAbstract.setBase - Index:'+ i2s(Index)+ '  NoOfBas:'+ i2s(NoOfBases));
 if CellsMX[Index].Base<> Value then
  begin
   CellsMX[Index].Base:= Value;
   DirtyBases:= TRUE;
  end;
end;



{--------------------------------------------------------------------------------------------------
   BASE EDIT

   To be used ONLY when the user manually replaces a base!
   It store the original base, replace it with the new base, then it redetects vectors.
--------------------------------------------------------------------------------------------------}
procedure TCubeAbstract.EditBase(CONST BasePos: integer; NewBase: TBase; CONST EditMxIndex: Integer);
begin
 { Store original base }
 if CellsMX[BasePos].EdtOrig= noBase
 then CellsMX[BasePos].EdtOrig:= Base[BasePos];

 { Store new base }
 Base[BasePos]:= NewBase;

 { Store edit point }
 CellsMX[BasePos].EditPtr:= EditMxIndex;                                                           {This makes sense only if this cube is a input sequence. If it is a contig, just set it to zero. }

 { Recompute vectors }
 DetectVectors;
end;


procedure TCubeAbstract.EditBase(CONST BasePos: integer; NewBase: TBase; CONST EditMxIndex: Integer; Mismatch: TMismStat);
begin
 { Store original base }
 if CellsMX[BasePos].EdtOrig= noBase
 then CellsMX[BasePos].EdtOrig:= Base[BasePos];

 { Store new base }
 Base[BasePos]:= NewBase;

 { Store edit point }
 CellsMX[BasePos].EditPtr:= EditMxIndex;                                                           {This makes sense only if this cube is a input sequence. If it is a contig, just set it to zero. }


 CellsMX[BasePos].MismStatus:= Mismatch;

 { Recompute vectors }
 DetectVectors;
end;




function TCubeAbstract.getEdtOrig(Index: integer): TBase;
begin
 Assert( (Index>= ctCellsIndex) AND (Index<= NoOfBases),'TCubeAbstract.getEdtOrig - Index:'+ i2s(Index)+ '  NoOfBas:'+ i2s(NoOfBases));
 Result:= CellsMX[Index].EdtOrig;
end;


procedure TCubeAbstract.setEdtOrig(Index: integer; Value: TBase);
begin
 Assert( (Index>= ctCellsIndex) AND (Index<= NoOfBases), 'TCubeAbstract.setEdtOrig - Index:'+ i2s(Index)+ '  NoOfBas:'+ i2s(NoOfBases));
 CellsMX[Index].EdtOrig:= Value;
end;




function TCubeAbstract.getEditPoint (Index: Integer): Integer;
begin
 Assert(Index>= ctCellsIndex, 'getEditPoint: ' + IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases   , 'getEditPoint_: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 Result:= CellsMX[index].EditPtr;
end;


procedure TCubeAbstract.setEditPoint (Index, Value: Integer);
begin
 Assert(Index>= ctCellsIndex, 'setEditPoint: ' + IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases   , 'setEditPoint_: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 CellsMX[index].EditPtr:= Value;
end;









{--------------------------------------------------------------------------------------------------
   SNPs
--------------------------------------------------------------------------------------------------}

function TCubeAbstract.CellIsIUPAC(Index: Integer): Boolean;
begin
 Result:= CharInSet(Base[Index], Ambiguity);
end;

function TCubeAbstract.CellIsIUPACDetected(Index: Integer): Boolean;
begin
 Result:= (Base[Index] <> Cell(Index).BaseOrig) AND CharInSet(Base[Index], Ambiguity);
end;







function TCubeAbstract.CellIsMutDiff (Index: integer): Boolean;                                    { Valid only for a contig base! Returns true if the contig base is an IUPAC base bacause the input sequences are different at this position }
begin
 Result:= CellsMX[Index].CtgMutation= cmDifference;
end;

function TCubeAbstract.CellIsMutIndel (Index: integer): Boolean;
begin
 Result:= CellsMX[Index].CtgMutation= cmIndel;
end;

function TCubeAbstract.CellIsMutant (Index: integer): Boolean;     { Returns TRUE if the cell has ANY kind of mutation (SNP, Indel, Difference. What about 'UserEdited' ????? It seems that I don't look at UserEdited  }
begin
 Result:= CellIsMutDiff(Index) OR CellIsMutIndel(Index) OR CellIsIUPAC(Index);
end;








function TCubeAbstract.CellSnpArea(Index: Integer): Byte;                { Mutant peak to reference peak area ratio. I need it for the Log }
begin
 Result:= CellsMX[Index].SnpAreaRatio;
end;


function TCubeAbstract.CellSnpOverlap(Index: Integer): Byte;
begin
 Result:= CellsMX[Index].SnpOverlap;
end;




{===============================================================================
   CELLS
===============================================================================}
function TCubeAbstract.Cell(CONST Index: integer): RCell5;
begin
 Assert(Index>= ctCellsIndex, 'Cell : '+ IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases   , 'Cell_: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 Result:= CellsMX[Index];
end;


procedure TCubeAbstract.CellSet (CONST Index: Integer; CONST NewValue: RCell5);
begin
 Assert(Index>= ctCellsIndex, 'Cellset : '+ IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases   , 'Cellset_: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 CellsMX[Index]:= NewValue;
 DirtyBases:= TRUE;
end;


function TCubeAbstract.CellProperties(CONST BasePos: Integer): string;
begin
 Result:= 'Base: '+ Base[BasePos]+ '.  Pos: '+ i2s(BasePos)+ '.  QV: '+ i2s(CellQV[BasePos]);
end;





{--------------------------------------------------------------------------------------------------
   CELLS - Bookmark, Mismatch, Color
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.getBookmark(Index: integer): Boolean;
begin
 Assert((Index>= ctCellsIndex) AND (Index<= NoOfBases), 'TCubeAbstract.getBookmark: '+ i2s(Index));
 Result:= CellsMX[index].Bookmark;
end;

procedure TCubeAbstract.setBookmark(Index: integer; Value: Boolean);
begin
 Assert((Index>= ctCellsIndex) AND (Index<= NoOfBases), 'TCubeAbstract.setBookmark: '+ i2s(Index));
 CellsMX[index].Bookmark:= Value;
end;


function  TCubeAbstract.getMismatch (Index: integer): TMismStat;
begin
 Assert((Index>= ctCellsIndex) AND (Index<= NoOfBases), 'TCubeAbstract.getMismatch: '+ i2s(Index)); { CellsMX E INDEXATA IN 1 }
 Result:= CellsMX[index].MismStatus;
end;

procedure TCubeAbstract.setMismatch (Index: integer; Value: TMismStat);
begin
 Assert((Index>= ctCellsIndex) AND (Index<= NoOfBases), 'TCubeAbstract.setMismatch: '+ i2s(Index));
 CellsMX[index].MismStatus:= Value;
end;


function  TCubeAbstract.getColorB (Index: integer): Tcolor;
begin
 Assert(Index>= ctCellsIndex, 'getColorB: ' + IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases   , 'getColorB : '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));                                                                        { CellsMX E INDEXATA IN 1 }
 Result:= CellsMX[index].ColorB;
end;

procedure TCubeAbstract.setColorB (Index: integer; Value: Tcolor);
begin
 Assert(Index>= ctCellsIndex, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 CellsMX[index].ColorB:= Value;
end;

function  TCubeAbstract.getColorT (Index: integer): Tcolor;
begin
 Assert(Index>= ctCellsIndex, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 Result:= CellsMX[index].ColorT;
end;

procedure TCubeAbstract.setColorT (Index: integer; Value: Tcolor);
begin
 Assert(Index>= ctCellsIndex, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 CellsMX[index].ColorT:= Value;
end;









{===============================================================================
   CELLS - CONVERSION - Base 2 Sample
===============================================================================}
function  TCubeAbstract.GetPtr2Smpl(Index: integer): Integer;                                      { Pointer to sample }  { Intoarce sample-ul corespunzator pt baza data ca parametru }
{$IFDEF CheckValidBase2Smpl}
VAR i, Rez: Integer; MsgError: string;
{$ENDIF CheckValidBase2Smpl}
begin
 Assert(Index<= NoOfBases   , 'GetPtr2Smpl: '+ IntToStr(Index)+ '/NoOfBases"'+ IntToStr(NoOfBases));
 Assert(Index>= ctCellsIndex, 'GetPtr2Smpl: '+ IntToStr(Index)+ '/CellsIndex:'+ IntToStr(ctCellsIndex));                                                                     { CellsMX E INDEXATA IN 1 }

 Result:= CellsMX[Index].Ptr2Smpl;
 Assert((Result>= 0) AND (Result<= NoOfSamples), 'GetPtr2Smpl: '+ IntToStr(Result)+ '/NoOfSamples: '+ IntToStr(NoOfSamples));

 { DEBUG }
 {$IFDEF CheckValidBase2Smpl}
 for i:= ctChromaIndex to NrOfSamples DO
  if Chroma[i].Ptr2Base= BaseIndex then
    begin
     Rez:= i;
     if Rez<> Result
     then Bip(800, 10);
     Exit;
    end;
 MsgError:= 'No sample assigned to this base!'
   +CRLF
   +CRLF+ 'DETAILS:'
   +CRLF+ '  File name: '+ ParentFileName
   +CRLF+ '  Base: '     + BASES[BaseIndex]+ ';  Position: '+ i2s(BaseIndex)
   +CRLF+ '  NoOfBases: '+ i2s(NoOfBases)
   +CRLF+ '  NrOfSamples: '+ i2s(NrOfSamples)
   +CRLF
   +CRLF+ 'This message appears usually when your chromatograms are invalid of incorrect formatted. Please report this issue to us and we will fix it in the next 48 hours.';
 MsgError:= ReplaceChar(MsgError, #0, '_');
 MesajErrDetail ParameterOrderChanged('TCubeAbstract.Base2Smpl', MsgError);
 {$ENDIF CheckValidBase2Smpl}
end;


procedure TCubeAbstract.SetPtr2Smpl(Index: integer; Value: Integer);
begin
 Assert(Index>= ctCellsIndex, 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(ctCellsIndex));
 Assert(Index<= NoOfBases   , 'xxxxxx: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfBases));
 CellsMX[index].Ptr2Smpl:= Value;
end;


procedure TCubeAbstract.SetReciprocPointers(BaseIndex, SampleIndex: Integer);
begin
  Chroma[SampleIndex].Ptr2Base:= BaseIndex;
  Base2Smpl[BaseIndex]:= SampleIndex;
end;





{--------------------------------------------------------------------------------------------------
   RESET CELL
--------------------------------------------------------------------------------------------------}
procedure RCell5.Clear;
begin
 Base        := '?';                                                                                 { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
 QV          := 0;
 MismStatus  := msNormal;                                                                            { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
 Trimmed     := FALSE;
 EdtOrig     := noBase;                                                                              { Base before user edit }
 EditPtr     := 0;                                                                                   { Stores the number of the edit that the user made in contig, edit which coresponds to this base (in imput sample) }    { Everytime I edit a cell I increment EditCounter and I assign it to CellMX.EditPoint. Thereofre, each edit in the Input Lines gets a unique (incrementing) number. The same for contig with the exception that when I edit a base in contig all bases above it get the same EditPoint }
 Bookmark    := FALSE;
 Vector      := FALSE;                                                                               { True if this cell has a Vector assigned to it }
 ColorB      := 0;                                                                                   { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
 ColorT      := 0;                                                                                   { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
 Ptr2Smpl    := NoneAssigned;                                                                        { position of the sample associated with this base }
 SnpOrig_    := noBase;                                                                              { If different than NoBase then it means that a SNP was detected here. It does not store the SNP base but the original base as imported from SCF/ABI file }
 SnpAreaRatio:= 0;                                                                                   { Mutant peak to reference peak area ratio (%). I need it for the RamLog }
 SnpOverlap  := 0;                                                                                   { Overlap between mutant peak and reference peak (%). I need it for the Log  }
 CellProb.A  := 0;                                                                                   { Probability of it being an A }
 CellProb.C  := 0;                                                                                   { Probability of it being an C }
 CellProb.G  := 0;                                                                                   { Probability of it being an G }
 CellProb.T  := 0;                                                                                   { Probability of it being an G }
end;



function RCell5.BaseOrig: TBase; { Show bases as in original file (before user editing, SNP detection, recomputing values for N bases, etc). This applies both the asembly grid and the chromatograms. }
VAR Max: Byte;
begin
 Max:= CellProb.A;
 Result:= 'A';

 if  (CellProb.C = Max)
 AND (CellProb.G = Max)
 AND (CellProb.T = Max)
 then Result:= 'N'
 else
  begin
   if CellProb.C > Max then
    begin
     Max:= CellProb.C;
     Result:= 'C';
    end;

   if CellProb.G > Max then
    begin
     Max:= CellProb.G;
     Result:= 'G';
    end;

   if CellProb.T > Max
   then Result:= 'T';
  end;
end;


procedure RCell5.Reverse;
VAR x: Integer;
begin
 Base:= BaseInverse(Base);

 { Swap probabilities }
 x := CellProb.A;
 CellProb.A:= CellProb.T;
 CellProb.T:= X;

 X := CellProb.C;
 CellProb.C:= CellProb.G;
 CellProb.G:= X;
end;





procedure TCubeAbstract.CellReset(CONST BasePos: integer);
begin
 Assert((BasePos>= ctCellsIndex) AND (BasePos<= NoOfBases), 'TCubeAbstract.CellReset: '+ i2s(BasePos));
 CellsMX[BasePos].Clear;
end;


procedure TCubeAbstract.CellMxReset;
VAR i: Integer;
begin
 for i:= ctCellsIndex to NoOfBases
  DO CellReset(i);
end;




{--------------------------------------------------------------------------------------------------
   QV
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.GoodBasesPercent: Real;                                                     { Cat % din baze sunt bune (GoodQV) dupa TrimEngine }
begin
 Assert(NoOfGoodBases > 0 , 'GoodBasesPercent: '+ IntToStr(NoOfGoodBases));
 Assert(NoOfBases     > 0 , 'GoodBasesPercent.NoOfBases'+ IntToStr(NoOfBases));
 Result:= ProcentRepresent(NoOfGoodBases, NoOfBases);
end;


function TCubeAbstract.NoOfGoodBases: integer;
begin
 Result:= LastGoodBase- GoodQVStart+ 1;
 Assert(Result> 0, 'TCubeAbstract.NoOfGoodBases <= 0');
 if Result<= 0
 then Result:= 1;                                                                                  { nu am voie sa am 0 baze bune ca pot sa fut proiectul cand incerc sa fac ceva de genul   X= NoOfBases/ NoOfGoodBases }
end;


function TCubeAbstract.HasAmbiguousBases: Boolean;
VAR
   CurBase: Integer;
   aBase: TBase;
begin
 Result:= FALSE;
 for CurBase:= ctCellsIndex to NoOfBases DO                                                        { Check the entire string for ambigous bases }
  begin
   aBase:= Base[CurBase];
   if CharInSet(aBase, Ambiguity) OR CharInSet(aBase, DNA_N)
   then EXIT(TRUE);
  end;
end;


function TCubeAbstract.Sample2QV(CONST Sample: integer): byte;                                     { intoarce QV-ul corespunzator pt sample-ul dat ca parametru }
begin
 Assert(HasChroma, 'Sample2QV.HasChroma: FALSE');
 Assert(Sample> 0, 'Sample2QV.Sample: '+ IntToStr(sample));
 Assert(Sample<= NoOfSamples, 'Sample2QV.Sample: '+ IntToStr(Sample)+ '/'+ IntToStr(NoOfSamples));
 Result:= CellQV[Chroma[sample].Ptr2Base];
end;



{-------------------------------------------------------------------------------             |BBBBBB|okokokokokokokok|BBBBBBB|
   GOOD QV START / END                                                                              |G1            G2|
--------------------------------------------------------------------------------                                                }
function TCubeAbstract.getGoodQVStart: Integer;
VAR CurBase: Integer;
begin
 Result:= -1;
 for CurBase:= ctCellsIndex to NoOfBases DO
   if NOT CellsMX[CurBase].Trimmed then
    begin
     Result:= CurBase;
     Break;
    end;

 Assert(Result < NoOfBases, 'All bases have been marked as ''low quality'' for sample '+ {ExtractOnlyName(FileNameParent)} ScreenName+ CRLF+ 'This will result in a null length sample!'+ CRLF+ 'getGoodQVStart');

 { Force at least 1 good base! }
 if Result= -1
 then Result:= NoOfBases DIV 2;  // Assert(Result>= ctCellsIndex, 'All bases have been marked as ''low quality'' for sample '+ ExtractOnlyName(FileNameParent)+ CRLF+ 'This will result in a null length sample!');
end;


procedure TCubeAbstract.setGoodQVStart(Value: Integer);
VAR CurBase: Integer;
begin
 Assert((Value>= ctCellsIndex) AND (Value< NoOfBases), 'Cannot mark all bases as ''low quality'' (sample '+ {ExtractOnlyName(FileNameParent)} ScreenName+')'+ CRLF+ 'This will result in a null length sample!');

 { Mark from the beginning of the cub to current position as untrusted }
 for CurBase:= ctCellsIndex to Value-1                                                             { -1 because I mark as grey all bases BEFORE the specified position }
  DO CellsMX[CurBase].Trimmed:= TRUE;

 { Mark from the current position to GoodQVEnd as trusted }
 for CurBase:= Value to LastGoodBase
   DO CellsMX[CurBase].Trimmed:= FALSE;

 DirtyGoodBases:= TRUE;
end;


function TCubeAbstract.getGoodQVEnd: Integer;                                                      { ASTA E LENTA! }
VAR CurBase: Integer;
begin
 Result:= -1;
 for CurBase:= NoOfBases downto ctCellsIndex DO
   if NOT CellsMX[CurBase].Trimmed then
    begin
     Result:= CurBase;
     Break;
    end;
 {TODO 3: Result:= GoodQVEnd}
 //Assert(Result<> -1, 'All bases have been found marked as ''low quality'' for sample '+ ExtractOnlyName(FileNameParent)+ CRLF+ 'This will result in a null length sample!');
 //Assert(Result >  1, 'All bases will be marked as ''low quality'' for sample '+ ExtractOnlyName(FileNameParent)+ CRLF+ 'This will result in a null length sample!'+ CRLF+ 'getGoodQVEnd');
end;


procedure TCubeAbstract.setLastGoodBase(LastGoodBase: Integer);                                    { Original name was GoodQVEnd }
VAR CurBase: Integer;
begin
 Assert((LastGoodBase>= ctCellsIndex) AND (LastGoodBase<= NoOfBases), 'Cannot mark all bases as ''low quality'' (sample '+ ScreenName {ExtractOnlyName(FileNameParent)}+')'+ CRLF+ 'This will result in a null length sample!');

 { Mark from GoodQVStart to current pos as trusted }
 for CurBase:= GoodQVStart to LastGoodBase
   DO CellsMX[CurBase].Trimmed:= FALSE;

 { Mark from the current position to end of string as untrusted }
 for CurBase:= LastGoodBase+1 to NoOfBases                                                         { +1 because I mark as grey all bases AFTER the specified position}
  DO CellsMX[CurBase].Trimmed:= TRUE;

 DirtyGoodBases:= TRUE;
 {TODO 3: GoodQVEnd:= LastGoodBase}
end;



{ GoodQVStart { including     |??????|kkkkkkkkkkkkkkkk|???????|                 ATENTIE!  E indexat in 1
  GoodQVEnd   { including           B|Good1      Good2|B                        }
function TCubeAbstract.getGoodStSmpl: Integer;                                                     { same as GoodQVStart but returns a sample nr instead of base nr }
begin
 Result:= Base2Smpl[GoodQVStart];
end;

function TCubeAbstract.getGoodEnSmpl: Integer;
begin
 Result:= Base2Smpl[LastGoodBase];
end;



function TCubeAbstract.getBaseGr(Index: integer): boolean;                                         { Returns true if the base is gray }
begin
 Assert(Index>= ctCellsIndex, 'TCubeAbstract.getBaseGr: '+ i2s(Index));
 Assert(Index<= NoOfBases,    'TCubeAbstract.getBaseGr: '+ i2s(Index));
 Result:= CellsMX[Index].Trimmed;
end;


procedure TCubeAbstract.setBaseGr(Index: Integer; Value: Boolean);
begin
 Assert(Index>= ctCellsIndex, 'TCubeAbstract.getBookmark: '+ i2s(Index));
 Assert(Index<= NoOfBases,    'TCubeAbstract.getBookmark: '+ i2s(Index));
 CellsMX[Index].Trimmed:= Value;
end;




{--------------------------------------------------------------------------------------------------
   CELL QV
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.getQV(Index: integer): Byte;
begin
 Assert(Index>= ctCellsIndex, 'TCubeAbstract.getBaseQV: '+ i2s(Index));
 Assert(Index<= NoOfBases   , 'TCubeAbstract.getBaseQV: '+ i2s(Index));
 Result:= CellsMX[Index].QV;                                                                       { CellsMX E INDEXATA IN 1 }
end;


procedure TCubeAbstract.setQV(Index: integer; Value: byte);                                        { Dupa asta NU trebuie sa apelez ApplyTrimEng, ca sa nu suprascriu editarile facute manual de user }
begin
 Assert((Index>= ctCellsIndex) AND (Index<= NoOfBases), 'TCubeAbstract.setBaseQV: '+ i2s(Index));
 if Value> 100 then Value:= 100;                                                                   { Check against wrong input (Tratez cazul Bradley 'PXY21-1_B_92294.scf' unde valorile lui QV poate sa fie pana la 255) }
 CellsMX[Index].QV:= Value;
end;




function TCubeAbstract.getCellProb(Index: integer): RBaseProb;
begin
 Assert(Index>= ctCellsIndex, 'TCubeAbstract.getCellProb: '+ i2s(Index));
 Assert(Index<= NoOfBases   , 'TCubeAbstract.getCellProb: '+ i2s(Index));
 Result:= CellsMX[Index].CellProb;                                                                 { CellsMX E INDEXATA IN 1 }
end;














{===============================================================================
   SAMPLES
===============================================================================}
procedure TCubeAbstract.setNrSamples(CONST Value: Integer);
CONST
   ctMaximumSequence= 30000000;                                                                    { Limit to 30 mil. samples. Echivalent of 1.5 MBases }
begin
 Assert(Value>= 0, 'setNrSamples.Value: '+ IntToStr(Value));
 Assert(Value<= ctMaximumSequence, 'The length of this sample exceeds the capacity of this program.');

 if FSamples<> value then
  begin
   FSamples:= Value;
   SetLength(Chroma, FSamples+ ctChromaIndex);
  end;
end;


function TCubeAbstract.getSample(Index: integer): RSample;
begin
 Assert(Index>= ctChromaIndex                , 'getSample.Index:  '+ IntToStr(Index)+ '/'+ IntToStr(ctChromaIndex));
 Assert(Index<= NoOfSamples-1 +ctChromaIndex, 'getSample.Index_: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfSamples));
 Result:= Chroma[index];
end;


procedure TCubeAbstract.setSample(Index: Integer; CONST Value: RSample);                           { Access to chromatogram's samples }
begin
 Assert(Index>= ctChromaIndex                 , 'setSample: '      + IntToStr(Index)+ '/'+ IntToStr(ctChromaIndex));
 Assert(Index<= (NoOfSamples-1) +ctChromaIndex, 'setSample.Index: '+ IntToStr(Index)+ '/'+ IntToStr(NoOfSamples));
 Chroma[index]:= Value;
end;




{--------------------------------------------------------------------------------------------------
   SAMPLES: CONVERSIONS
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.Sample2Base(CONST Sample: integer): TBase;                                  { intoarce baza corespunzatoare pt sample-ul dat ca parametru }  {recent added}
begin
 Assert(HasChroma, 'Sample2Base.HasChroma: false');
 Assert(Sample>= ctChromaIndex, 'Sample2Base.Sample: '+ IntToStr(Sample));
 Assert(Sample<= NoOfSamples  , 'Sample2Base.Sample: '+ IntToStr(Sample)+ '/'+ IntToStr(NoOfSamples));

 if SampleHasBaseAssigned(Sample)
 then Result:= Base[Chroma[Sample].Ptr2Base]
 else Result:= noBase;
end;


function TCubeAbstract.Sample2BaseI(CONST Sample: integer): integer;                               { Base pos corresponding to this sample pos }
begin
 Assert(HasChroma);
 Assert(Sample> 0           , 'TCubeAbstract.Sample2BaseI: '+ i2s(Sample));
 Assert(Sample<= NoOfSamples, 'TCubeAbstract.Sample2BaseI: '+ i2s(Sample)+ '/'+ IntToStr(NoOfSamples));

 Result:= Chroma[Sample].Ptr2Base;
end;


function TCubeAbstract.Sample2BaseApprox(CONST Sample: integer): integer;                          { Base pos corresponding to this sample pos. If there is no base assigned to this sample, than find the closes neigbor }
VAR i, Left, Right: Integer;
begin
 Assert(Sample> 0);

 if Sample> NoOfSamples                                                                            { Asta era un Assert, dar asa e mai bine }
 then EXIT(NoOfBases);

 Result:= -1;

 if Chroma[Sample].Ptr2Base<> NoneAssigned
 then Result:= Chroma[Sample].Ptr2Base
 else
  begin
   Left := Sample;
   Right:= Sample;
   for i:= 0 to High(Integer) DO                                                                   { Initially, both the start at the center. The 'cursor' is then moved apart from the center until a base is found }
    begin
     Dec(Left);
     Inc(Right);
     if (Right> NoOfSamples) OR (Left< ctCellsIndex)
     then EXIT(-1);

     if (Chroma[Left].Ptr2Base<> NoneAssigned)
     then EXIT(Chroma[Left].Ptr2Base);
     if (Chroma[Right].Ptr2Base<> NoneAssigned)
     then EXIT(Chroma[Right].Ptr2Base);
    end;
  end;
end;


function TCubeAbstract.Base2SampleHeight(CONST BaseNo: Integer): Integer;                          { Gets a base as input and returns the sample height for the peak associated to that base }
begin
 Result:= Base2SampleHeight(Base[BaseNo], Base2Smpl[BaseNo]);
end;


function TCubeAbstract.Base2SampleHeight(CONST BaseNo, aSample: Integer): Integer;                 { Gets a base position as input and returns the sample height for the peak associated to that base }
begin
 Result:= Base2SampleHeight(Base[BaseNo], aSample);
end;


function TCubeAbstract.Base2SampleHeight(CONST aBase: TBase; aSample: Integer): Integer;           { Gets a base as input and returns the sample height for the peak associated to that base }
begin
 case UpCase(aBase) of
   'A': Result:= Sample[aSample].HeightA;
   'C': Result:= Sample[aSample].HeightC;
   'G': Result:= Sample[aSample].HeightG;
   'T': Result:= Sample[aSample].HeightT;
  else
    Result:= 0;                                                                                    {BASE DOES NOT HAVE TRACE! This happens for N, and IUPAC bases! }
 end;
end;




{--------------------------------------------------------------------------------------------------
   SAMPLES - CONVERSIONS
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.SampleHasBaseAssigned(CONST Sample: integer): Boolean;                      { This checks if the specified sample has a base assigned to it (doesn't matter which color is that base) }
begin
 Result:= Chroma[Sample].Ptr2Base> 0;
end;


function TCubeAbstract.SampleHasBaseAssigned(CONST Sample: Integer; CONST BaseColor: TBase): Boolean;  { This checks if the specified sample has a base of the specified color assigned to it }
begin
 Result:= SampleHasBaseAssigned(Sample);
 if Result
 then Result:= Base[Chroma[Sample].Ptr2Base]= BaseColor;
end;


function TCubeAbstract.SampleHasBaseAssignedEx(CONST Sample: Integer; CONST BaseColor: TBase): Boolean;  { Same as above but it ignores the BaseColor parameter to Ns. In other words it always returns TRUE if there is a base at this sample and that base is N. }
VAR CurBase: TBase;
begin
 Result:= SampleHasBaseAssigned(Sample);

 if Result then
  begin
   CurBase:= Sample2Base(Sample);
   if CharInSet(CurBase, DNA_N)                                                                             { Ignore the BaseColor parameter if the assigned base is N }
   then Result:= TRUE
   else Result:= Base[Chroma[Sample].Ptr2Base]= BaseColor;
  end;
end;


function TCubeAbstract.SampleHasBaseAssignedArround(CONST Sample: Integer): Boolean;               { I check few positions arround this sample because for some strange reasons, in ABI, the peaks are shifted to the left with 1 or 2 pixels }
begin
 Assert(HasChroma);
 Result:= (
   (Chroma[Sample+1].Ptr2Base> 0) OR
   (Chroma[Sample  ].Ptr2Base> 0) OR
   (Chroma[Sample-1].Ptr2Base> 0) OR
   (Chroma[Sample-2].Ptr2Base> 0) OR
   (Chroma[Sample-3].Ptr2Base> 0) );
end;


procedure TCubeAbstract.Ulimination;                                                               { Replaces the 'U' bases with 'T' }
VAR s: BaseString;
begin
 s:= bases;
 ReplaceChar(s, 'U', 'T');
 ReplaceChar(s, 'u', 'T');
 setBases(s);
end;





{--------------------------------------------------------------------------------------------------
   CHROMA
--------------------------------------------------------------------------------------------------}
function TCubeAbstract.HasChroma;
begin
 Result:= Length(Chroma) > 0;
end;


function TCubeAbstract.HighestTrace(aSample: Integer): Integer;                                    { Check the ACGT traces for the highest value }
begin
 Result:= 0;
 if Sample[aSample].HeightA> Result then Result:= Sample[aSample].HeightA;
 if Sample[aSample].HeightC> Result then Result:= Sample[aSample].HeightC;
 if Sample[aSample].HeightG> Result then Result:= Sample[aSample].HeightG;
 if Sample[aSample].HeightT> Result then Result:= Sample[aSample].HeightT;
end;


function TCubeAbstract.DistanceBetween(CONST Base1, Base2: Integer): Integer;                      { The distance (in samples) between two bases }
begin
 Assert(Base1< Base2);
 Assert(Base1>= ctCellsIndex, i2s(Base1));
 Assert(Base2<= NoOfBases );

 Result:= Base2Smpl[Base2] - Base2Smpl[Base1];
end;


procedure TCubeAbstract.ExtractTrace(TraceColor: TBase; OUT Trace: RTrace);                        { Returns the content of the specified trace (ACGT) }
VAR smpl: Integer;
begin
 TraceColor:= UpCase(TraceColor);
 Assert(CharInSet(TraceColor, DNA_ACGT));                                                                   { This is because ExtractTrace can only work on one of the 4 traces. 'N' bases are not connected to any of the 4 traces }

 SetLength(Trace.Height, Length(Chroma));                                                          { Make it of the SAME length }
 Trace.Color:= TraceColor;

 for smpl:= 0 to NoOfSamples-ctChromaIndex DO                                                      { Copy entierly even if the first element has no meaningful value in it }
   case TraceColor of
     'A': Trace.Height[smpl]:= Chroma[smpl].HeightA;
     'C': Trace.Height[smpl]:= Chroma[smpl].HeightC;
     'G': Trace.Height[smpl]:= Chroma[smpl].HeightG;
     'T': Trace.Height[smpl]:= Chroma[smpl].HeightT;
   end;
end;


function TCubeAbstract.AveragePeakHeight: Integer;                                                 { Media peak-urilor (care au baze asignate) pt toata cromatograma }
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











{===============================================================================
   SOFT EDIT - INSERT BASE
================================================================================
 Description:
     Insereaza intre doua baze deja existete
     BeforeBase - arata baza din UnitsMX DUPA care baza va fi inserata baza curenta
     Nu pot sa introduc prea multe baze una dupa alta (limita e cam 6) pt ca
       devin prea inghesuite (spatiul dintre ele = 1 pixel, si nu mai am
       unde sa o introduc si pe a 7-a
===============================================================================}
procedure TCubeAbstract.InsertBase (CONST NewBase: TBase; CONST QV: Byte; CONST BeforeBase: Integer);  { NEW } { The base at position 'BeforeBase' and all bases after this position will be pushed to the right. The empty placed created at 'BeforeBase' position will be filled by the new base }
VAR
   MustCreate: Boolean;
   i, InsertAtSample, AfterBase: Integer;
begin
 if BeforeBase= 1 then
  begin
   AppendBaseAtBeg(NewBase, QV);                                                                   { Daca vreau sa introduc o baza noua pe prima pozitie, atunci in stanga ei nu am nici un vecin }
   EXIT;
  end;

 InsertAtSample:= -777;                                                                            { Shut up the compiler }
 AfterBase:= BeforeBase-1;

 { Add one extra base }
 NoOfBases:= NoOfBases+ 1;                                                                         { rezerv memorie de pe acum pt ca am nevoie de ea cand apelez Changed sa refac bazele }

 if HasChroma then
  begin
   { Where can I insert? }
   InsertAtSample:= GetMiddlePoint(AfterBase, MustCreate);                                         { Returns the midle distance between two bases. If MustCreate=true I need to insert a new sample right after the value returned. }

   if MustCreate then
    begin
     Inc(InsertAtSample);                                                                          { If MustCreate=true I need to insert a new sample right after the value returned. }
     NoOfSamples:= NoOfSamples+ 1;                                                                 { Enlarge matrix }
     { Shift all samples }
     for i:= NoOfSamples-1 downto InsertAtSample                                // because I have +1 below
       DO Chroma[i+1]:= Chroma[i];
    end;

   { Shift all pointers }
   for i:= InsertAtSample to NoOfSamples DO
     if Chroma[i].Ptr2Base> NoneAssigned
     then Inc(Chroma[i].Ptr2Base);

   { Setup new sample }
   Chroma[InsertAtSample].Ptr2Base:= BeforeBase;
  end;

 { Shift all cells }
 for i:= NoOfBases-1 downto BeforeBase
  DO CellSet(i+1, Cell(i));

 { Shift pointers from bases to peaks }                                         { needed because I entered a new sample and now all pointers are behind with 1 }
 if MustCreate then
  for i:= NoOfBases-1 downto BeforeBase
   DO CellsMX[i].Ptr2Smpl:= CellsMX[i].Ptr2Smpl + 1;

 { Insert the new base }
 CellReset(BeforeBase);                                                                            { I need this in order to initialize the cell. Among others it will set CellEditPoint:= 0 }
 Base     [BeforeBase]:= NewBase;
 CellQV   [BeforeBase]:= QV;
 CellTrimmed[BeforeBase]:= CellTrimmed[BeforeBase-1];                                              { Copy value from its neighbor. This is very important else I will enter a trimmed (gray) base in the middle of the GoodBases segment (or vice versa) }

 if HasChroma
 then SetReciprocPointers(BeforeBase, InsertAtSample);

 DirtyBases:= TRUE;
end;



function TCubeAbstract.GetMiddlePoint (CONST FirstBase: Integer; OUT MustCreate: Boolean): Integer;  { Returns the midle distance between two bases. If MustCreate=true I need to insert a new sample right after the value returned. }
VAR SmplVecinSt, SmplVecinDr: Integer;
begin
 Assert(HasChroma);                                                                                { I don't have to call this function because I CAN ALLWAYS INSERT BASES IN FASTA/TXT FILES DIRECTLY }
 Assert(FirstBase> 0, i2s(FirstBase));
 Assert(FirstBase<= NoOfBases, i2s(FirstBase));

 { Vecin stanga }
 SmplVecinSt:= Base2Smpl[FirstBase];

 { Vecin dreapta }
 if FirstBase+ 1 > NoOfBases                                                                       { Afla unde e prima baza in dreapta }
 then SmplVecinDr:= NoOfSamples                                                                    { I don't have anymore bases to the right }
 else SmplVecinDr:= Base2Smpl[FirstBase+1];

 { Compute middle point }
 Result:= SmplVecinSt+ (SmplVecinDr - SmplVecinSt) DIV 2;
 MustCreate:= (Result= SmplVecinSt) OR (Result= SmplVecinDr);
end;













{===============================================================================
   SOFT EDIT - APPEND BASE
===============================================================================}
procedure TCubeAbstract.AppendBaseAtBeg (CONST NewBase: TBase; CONST QV: Byte);                    { Append base at the beginning of the chromatogram }
VAR i: Integer;
begin
 if HasChroma then
  begin
    { CHROMA - Shift all samples }
    NoOfSamples:= NoOfSamples+ 12;
    for i:= ctChromaIndex To NoOfSamples DO
      if Chroma[i].Ptr2Base> NoneAssigned
      then Inc(Chroma[i].Ptr2Base);
   end;

 { BASES - Shift all cells }
 NoOfBases:= NoOfBases+ 1;
 for i:= NoOfBases-1 DownTO 1                                                                      { from NoOfBases to the Insertion position }
   DO CellSet(i+1, Cell(i));

 { Insert the new base }
 CellReset(1);
 Base[1]:= NewBase;
 CellQV[1]:= QV;
 CellTrimmed[1]:= CellTrimmed[2];                                                                  { Copy value from its neighbor.          <- this fixes: When I enter a new base in a ABI file right on position 1, I get an Assertion failure in BuildBases }

 if HasChroma
 then SetReciprocPointers(1, 1);

 DirtyBases:= TRUE;
 DetectVectors;
end;


procedure TCubeAbstract.AppendBaseAtEnd2 (CONST NewBase: TBase; CONST QV: Byte);                   { Append base at the end of the chromatogram }  { Important: The caller should call DetectVectors after that }
VAR i, InsertAtSample: Integer;
begin
 InsertAtSample:= -777;                                                                            { To make the compiler shutup }

 { Expand cub }
 NoOfBases:= NoOfBases+ 1;

 { CHROMA }
 if HasChroma then
  begin
   NoOfSamples:= NoOfSamples+ 24;
   InsertAtSample:= NoOfSamples- 12;
   Chroma[InsertAtSample].Ptr2Base:= NoOfBases;
   for i:= NoOfSamples- 24 to NoOfSamples
    DO Chroma[i].Ptr2Base:= NoneAssigned;                                                          { Initialize pointers }
  end;

 { BASES }
 CellReset(NoOfBases);
 Base     [NoOfBases]:= NewBase;
 CellQV   [NoOfBases]:= QV;
 CellTrimmed[NoOfBases]:= CellTrimmed[NoOfBases-1];                                                { Copy value from its neighbor }

 if HasChroma
 then SetReciprocPointers(NoOfBases, InsertAtSample);

 DirtyBases:= TRUE;
end;


procedure TCubeAbstract.AppendCellAtEnd (CONST NewCell: RCell5);                                    { Append base at the end of the chromatogram }  { Important: The caller should call DetectVectors after that }
VAR i, InsertAtSample: Integer;
begin
 InsertAtSample:= -777;                                                                            { To make the compiler shutup }

 { Expand cub }
 NoOfBases:= NoOfBases+ 1;

 { CHROMA }
 if HasChroma then
  begin
   NoOfSamples:= NoOfSamples+ 24;
   InsertAtSample:= NoOfSamples- 12;
   Chroma[InsertAtSample].Ptr2Base:= NoOfBases;
   for i:= NoOfSamples- 24 to NoOfSamples
    DO Chroma[i].Ptr2Base:= NoneAssigned;                                                          { Initialize pointers }
  end;

 { BASES }
 CellSet(NoOfBases, NewCell);
 CellTrimmed[NoOfBases]:= CellTrimmed[NoOfBases-1];                                                { Copy value from its neighbor }

 if HasChroma
 then SetReciprocPointers(NoOfBases, InsertAtSample);

 DirtyBases:= TRUE;
end;






{===============================================================================
   SOFT EDIT - REPLACE BASE
===============================================================================}
procedure TCubeAbstract.ReplaceBase;                                                               { Replaces a base in CellsMX si BASES. Dupa asta NU trebuie sa apelez ApplyTrimEng, ca sa nu suprascriu editarile facute manual de user }
begin
 if CellQV [BasePos]<> QV
 then CellQV [BasePos]:= QV;                                                                       { replace QV }

 if Base[BasePos] <> NewBase
 then Base[BasePos]:= NewBase;                                                                     { This will call "Changed" and it will set 'CellsMX[Position].Edited' to 'TRUE' }
end;


procedure TCubeAbstract.ReplaceBase2Gap(CONST BasePos: integer);                                   { Replace the specified base with a gap }
begin
 CellQV[BasePos]:= 0;
 Base  [BasePos]:= GAP;                                                                            { This will call "Changed" }
end;





{===============================================================================
   SOFT EDIT - OTHER EDITS
===============================================================================}
procedure TCubeAbstract.ReinsertGoodBases(CONST NewMiddle: BaseString);                            { REDISTRIBUI QV-urile  -  trebuie sa reasociez QV-urile cu bazele din cauza ca acum sirul de baze s-a modificat (a devenit mai lung prin adaugarea de noi GAP-uri) }
VAR
   i, Col: Integer;
   Rebuilt, LowQV1, LowQV2: BaseString;
begin
 { Detect gray ends }
 for i:= 1 to GoodQVStart-1
   DO LowQV1:= LowQV1+ Base[i];                                                                    { Capetele GRI nu au intrat in asamblare asa ca trebuie sa le determin }

 for i:= LastGoodBase+1 to NoOfBases
   DO LowQV2:= LowQV2+ Base[i];

 Rebuilt:= LowQV1+ NewMiddle+ LowQV2;

 { Replace GapAsm with real Gaps }
 for col:= 1 to Length(Rebuilt) DO
  if Rebuilt[col]= GapAsm
  then InsertBase(Gap, 0, col);                                                                    {! This will call also "DirtyBases" }  { The base at position 'BeforeBase' and all bases after this position will be pushed to the right. The empty placed created at 'BeforeBase' position will be filled by the new base }
end;


procedure TCubeAbstract.SwapCells(CONST Cell1, Cell2: Integer);                                    { UNUSED FUNCTION }
VAR Cel: RCell5;
begin
 Assert((Cell1> 0) AND (Cell1<= NoOfBases));
 Assert((Cell2> 0) AND (Cell2<= NoOfBases));

 Cel:= Cell(Cell1);
 CellSet(Cell1, Cell(Cell2));
 CellSet(Cell2, Cel);
 { DirtyBases:= TRUE --> this is called anyway by CellSet }
end;





{===============================================================================
   SOFT EDIT - DELETE BASES
===============================================================================}

{ Delete bases from the beginning of the chromatogram to BasePos (including).
  The resulted chromatogram will be shorter }
procedure TCubeAbstract.TruncateTo (CONST BasePos: Integer);
VAR i, AtSample: integer;
begin
 AtSample:= 0;

 if HasChroma then
   begin
    AtSample:= Base2Smpl[BasePos] +1;   { +1 in order to cut also the sample assigned to Base[BasePos] that will be deleted }

    { Reindex pointers to bases }
    for i:= AtSample {ChromaIndex} to NoOfSamples DO
     begin
      if   Chroma[i].Ptr2Base> NoneAssigned
      then Chroma[i].Ptr2Base:= Chroma[i].Ptr2Base- BasePos;
     end;

    { Shift samples }
    for i:= ChromaIndex to NoOfSamples- AtSample DO
     begin
       Chroma[i]:= Chroma[i+ AtSample];
     end;

    { Resize chromatogram }
    NoOfSamples:= NoOfSamples- AtSample;
   end;

 { Shift bases x positions to the left }
 for i:= CellsIndex to NoOfBases-BasePos
  DO CellSet(i, Cell(i+BasePos));                                { mut toate bazele cu o pozitie mai in fata }

 { Shift pointers to samples }
 if HasChroma then
  for i:= CellsIndex to NoOfBases-basepos DO
   begin
    if CellsMX[i].Ptr2Smpl > NoneAssigned
    then CellsMX[i].Ptr2Smpl:= CellsMX[i].Ptr2Smpl- AtSample;
   end;

 { Resize CellMX }
 NoOfBases:= NoOfBases-BasePos;
 DirtyBases:= TRUE;
end;



{ Delete bases from BasePos (including) to the end of the chromatogram.
  The resulted chromatogram will be shorter }
procedure TCubeAbstract.TruncateFrom (CONST BasePos: Integer);
VAR i, AtSample, LeftBases: integer;
begin
 LeftBases:= BasePos-1;

 if HasChroma then
   begin
    AtSample:= Base2Smpl[BasePos] -1;  { -1 in order to cut also the sample assigned to Base[BasePos] that will be deleted }

    { Resize chromatogram }
    NoOfSamples:= AtSample;

    { Check pointers }
    for i:= ChromaIndex to NoOfSamples DO
     begin
      if Chroma[i].Ptr2Base>= BasePos
      then RAISE Exception.Create('Bad chroma end!') { Chroma[i].Ptr2Base:= LeftBases }
     end;
   end;

 { Resize CellMX }
 NoOfBases:= LeftBases;

 { Shift pointers to samples }
 if HasChroma then
  for i:= CellsIndex to NoOfBases DO
   begin
    if CellsMX[i].Ptr2Smpl > NoOfSamples
    then RAISE Exception.Create('Bad CellsMX end!');   { CellsMX[i].Ptr2Smpl:= NoOfSamples }
   end;

 DirtyBases:= TRUE;
end;



{ Removes a base from the midle of a seq.
  The peak assigned to that base will remain there (the chromatogram will not be shortened) }
procedure TCubeAbstract.RemoveBase (CONST BasePos: Integer);
VAR i, AtSample: integer;
begin
 if HasChroma then
   begin
    { Sterg baza din chromatograma (Chroma) }
    AtSample:= Base2Smpl[BasePos];
    Chroma[AtSample].Ptr2Base:= NoneAssigned;                                                      { Sterge Baza }

    { Reindexeaza pointerii }
    for i:= AtSample+1 to NoOfSamples DO
      if   Chroma[i].Ptr2Base> NoneAssigned
      then Dec(Chroma[i].Ptr2Base);
   end;

 { Shift pointers }
 { Not needed because the samples (peak) remains there }

 { Shiftez toate bazele in UnitsMX one pos to the left }
 for i:= BasePos to NoOfBases-1
  DO CellSet(i, Cell(i+1));                                                                        { mut toate bazele cu o pozitie mai in fata }

 { Redimensionez CellMX }
 NoOfBases:= NoOfBases-1;                                                                          { asta o sa faca SetLength(UnitsMX, Numar_of_Bases) }

 DirtyBases:= TRUE;
end;


function TCubeAbstract.RemoveBasePossible (CONST BasePos: Integer): Boolean;                       { is it possible to delete this base? }
begin
 Result:= (NoOfBases>= ctMinimimSequence)
      AND (BasePos<= NoOfBases)
      AND (BasePos> 0);                                                                            { don't delete any more bases if I have only 10 bases left }
end;



procedure TCubeAbstract.RemoveAllGaps;                                                             {TODO 5: improve the algorithm to remove all bases at once }
VAR CurBase: Integer;
begin
 for CurBase:= NoOfBases-1 downto 1 DO
  if (Base[CurBase]= gap)
  then RemoveBase(CurBase);
end;














{===============================================================================
   REVERSE COMPLEMENT
===============================================================================}
procedure TCubeAbstract.Reverse;
VAR
  CurBase, CurSmpl, BzNr: Integer;
  Smpl       : RSample;
  tmpChroma  : AChromaMX;                                                                          { indexata in 1 }
  tmpCellMX  : ACells;
  ValEsantion: Word;
begin
 FReversed:= NOT FReversed;

 (* Reverse the chromatogram *)
 if HasChroma then
  begin
   BzNr:= -1;
   SetLength(tmpChroma, Length(Chroma));                                                           { Make it of the SAME length }
   for CurSmpl:= chromaindex to NoOfSamples DO                                                               { IGNORE the first element because it has no meaningful value in it }
    begin
     { extrag un sample si ii inversez matele }
     Smpl:= Chroma[CurSmpl];
     if Smpl.Ptr2Base> NoneAssigned then
      begin
       inc(BzNr);
       Smpl.Ptr2Base:= NoOfBases- BzNr;
      end;

     { inversez cele 4 sample-uri intre ele }
     ValEsantion := Smpl.HeightA;
     Smpl.HeightA:= Smpl.HeightT;
     Smpl.HeightT:= ValEsantion;

     ValEsantion := Smpl.HeightC;
     Smpl.HeightC:= Smpl.HeightG;
     Smpl.HeightG:= ValEsantion;

     { pun la loc in Chroma sampleurile, dar in ordine inversa }
     tmpChroma[NoOfSamples- CurSmpl+ 1]:= Smpl;
    end;
   Chroma:= tmpChroma;
  end;

 (* INTOARECE BAZELE *)
 BzNr:= -1;
 SetLength(tmpCellMX, Length(CellsMX));                   { I store the cells in reverse order in a temp matrix }
 for CurBase:= ctCellsIndex to NoOfBases DO
  begin
   inc(BzNr);
   tmpCellMX[CurBase]:= Cell(NoOfBases- BzNr);
   tmpCellMX[CurBase].Reverse;

   if HasChroma
   then tmpCellMX[CurBase].Ptr2Smpl:= NoOfSamples- tmpCellMX[CurBase].Ptr2Smpl+ 1;
  end;
 CellsMX:= tmpCellMX;           { Copy temp matrix back to the CellMX }

 DirtyBases:= TRUE;
 buildbases;                                                                                       { Force rebuild bases and goodbases. Fara asta NU merge. Nu stiu de ce }
 detectvectors;
end;


function TCubeAbstract.GoodBasesRC: BaseString;
begin
 Result:= ReverseComplement(GoodBases);
end;







{==================================================================================================
   VECTORS

   NOTE: Vectors are detected everytime the cub is edited/changed by calling 'DetectVectors' in 'getBases' function.
   WEIRD STUFF TO REMEMBER: If the vector contains an IUPAC base it will work (the vector will be detected and cut) because we extrapolate the vectors. However, if the CONTIG contains an IUPAC base in the vector region, the vector won't be detected because... well... there is no match between contig and the vector. That IUPAC base breaks the match.
==================================================================================================}
procedure TCubeAbstract.DetectVectors;
VAR cl: Integer;
begin
 inherited;   { <----- This is the real thing. It calls TVectorDetector.DetectVectors }

 { Mark bases in blue for left vector }
 if LeftVectorExists then
  for cl:= Vectors.Left.Starts TO Vectors.Left.Ends
    DO CellHasVector[cl]:= TRUE;

 { Mark bases in blue for right vector }
 if RightVectorExists then
  for cl:= Vectors.Right.Starts TO Vectors.Right.Ends
    DO CellHasVector[cl]:= TRUE;
end;


procedure TCubeAbstract.VectorsClearColor;
VAR BasePos: Integer;
begin
 inherited;

 if NoOfBases> 0 then                                                                              { An uninitialized object may have zero bases }
  for BasePos:= ctCellsIndex to NoOfBases                                                          { All existent cells already marked as Vectors, are cleared. }
   DO CellHasVector[BasePos]:= FALSE;
end;


function TCubeAbstract.GoodBasesNoVector: BaseString;
VAR iFrom, iTo: Integer;
begin
 CleanFromTo(iFrom, iTo);
 Result:= CopyTo(Bases, iFrom, iTo);
end;


function TCubeAbstract.VectorAtBase(BasePos: Integer): string;                                     { Returns the name of the vector associated with the specified base }
begin
 if (BasePos>= Vectors.Left.Starts) AND (BasePos<= Vectors.Left.Ends)
 then Result:= Vectors.Left.Name
 else
   if (BasePos>= Vectors.Right.Starts) AND (BasePos<= Vectors.Right.Ends)
   then Result:= Vectors.Right.Name;
end;


function TCubeAbstract.LeftVectorExists: Boolean;
begin
 Result:= Vectors.Left.Starts> 0;
end;


function TCubeAbstract.RightVectorExists: Boolean;
begin
 Result:= Vectors.Right.Starts> 0;
end;


function TCubeAbstract.BaseIsStriked(BasePos: Integer): Boolean;                                   { Returns true if the specified base will be cut from contig because it is located in vector or previous to the vector }
begin
 if TVectorDetector(Vectors.Detector).Keep
 then Result:= (LeftVectorExists AND (BasePos <  Vectors.Left.Starts)) OR (RightVectorExists AND (BasePos>  Vectors.Right.Ends))
 else Result:= (LeftVectorExists AND (BasePos <= Vectors.Left.Ends))   OR (RightVectorExists AND (BasePos>= Vectors.Right.Starts));
end;


function TCubeAbstract.getHasVect(Index: integer): Boolean;
begin
 Assert(Index>= ctCellsIndex);
 Assert(Index<= NoOfBases);
 Result:= CellsMX[index].Vector;
end;


procedure TCubeAbstract.setHasVect(Index: integer; Value: Boolean);
begin
 Assert(Index>= ctCellsIndex);
 Assert(Index<= NoOfBases);
 CellsMX[index].Vector:= Value;
end;


procedure TCubeAbstract.CleanFromTo(OUT FromBase, ToBase: Integer);                                      { Shows from which to which base I have to cut in order to remove the vectors and the low quality ends }
VAR Keep: Boolean;
begin
{ How it works:
    No vector found: Clean based only on GoodBases.
       Vector found: Clean based on GoodBases and Vector (whichever is found first)  }
 if Vectors.Detector = NIL
 then Keep:= TRUE
 else Keep:= TVectorDetector(Vectors.Detector).Keep;


 if Keep
 then
  begin
    { Left side of the sample }
    if (Vectors.Left.Starts< 1)                                                                    { Left vector does not exist (has not been detected) }
    OR (Vectors.Left.Starts< GoodQVStart)                                                          { If the vector starts/ends before the bad ends, then I ignore vector positon and copy based on quality of the ends. The thing is to throw away low quality ends. }
    then FromBase:= GoodQVStart
    else FromBase:= Vectors.Left.Starts;

    { Right side of the sample }
    if (Vectors.Right.Ends< 1)                                                                     { Right vector does not exist (has not been detected) }
    OR (Vectors.Right.Ends> LastGoodBase)
    then ToBase:= LastGoodBase
    else ToBase:= Vectors.Right.Ends;
  end
 else
  begin
    { Left side of the sample }
    if (Vectors.Left.Ends< 1)                                                                      { Left vector does not exist (has not been detected) }
    OR (Vectors.Left.Ends< GoodQVStart)
    then FromBase:= GoodQVStart
    else FromBase:= Vectors.Left.ends+1;                                                              { +1 because I want to copy bases AFTER the last base that is a vector }

    { Right side of the sample }
    if (Vectors.Right.Starts< 1)                                                                   { Right vector does not exist (has not been detected) }
    OR (Vectors.Right.Starts> LastGoodBase)
    then ToBase:= LastGoodBase
    else ToBase:= Vectors.Right.Starts-1;                                                             { -1 because I want to copy bases BEFORE the last base that is a vector }
  end;
end;


function TCubeAbstract.Search(CONST SearchStr: BaseString; Offset: Integer): integer;
begin
 Result:= PosEx(SearchStr, UpperCase(Bases), Offset);                                              { I need tu use UpperCase because some samples may contain low case bases }
 {TODO 5: Search: IMPROVE SPEED HERE! }
end;





END.

{ Documentation about CONSTRUCTORS: http://stackoverflow.com/questions/772336/using-inherited-in-the-create-constructor-of-an-tobject }



