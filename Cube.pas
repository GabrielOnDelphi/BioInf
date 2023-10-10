
UNIT Cube;

{==================================================================================================
 Heracle BioSoft SRL
 2016.11.01

 Object capabilities:
    + AsmGrid capabilities (Contig, Neighbours, ContigQuality)
    + Visual  capabilities (Colors)
    + QV graph

 Object tree:
    * CubeBase -> CubeBaseQV -> CubeBaseSNP -> CubeImporter -> Cube -> CubeStream
   
==================================================================================================}

INTERFACE
USES
   Winapi.Windows, System.SysUtils, System.Classes, Vcl.Graphics,
   ccCore, ccINIFile, CubeImporter, {cGraphics,} CubicDNA, clRamLog, ccRichLog;

TYPE
 TCubeObj = class(TCubeImport)
  private
    FColorBkg        : TColor;                                                                    { Color for normal cells (cells that are neither vectors (blue), erros (red), etc. This will give the overal color of the grid }
    FSolvedBkg       : TColor;
    FAutoContigColor : TColor;
    FColorESolvedGray: TColor;                                                                    { The color for solved background when it overlaps with Gray cell }
    FRainbowText     : TColor;                                                                    { Color of the text in rainbow mode }
    FRainbowBkg      : Boolean;
    FAsmOffset       : Integer;                                                                   { AsmGridOffset = Offset in AsmGrid. Indexed in 0. The TContigGrid will have to add 1 to this value, because on the position 0, it has the header. The AsmGrid don't have to add 1 because it is also indexed in 0. }
    procedure setColorBkg(const Value: TColor);
    procedure setSolvedBkg(const Value: TColor);
  protected
    function  getFileName: string; override;
    function  getBaseColor(Base: TBase): TColor;
    procedure setAsmOffset(Value: Integer);
  public
    {VISUAL/GRID RELATED}
    ColorEErrorB     : TColor;
    ColorBookmarkB   : TColor;
    FontColor        : TColor;                                                                    { Font color for any text (other than ACGTN-) }
    FBookmarkT       : TColor;
    BaseColorC       : Tcolor;                                                                    { bases }
    BaseColorG       : Tcolor;                                                                    { bases }
    BaseColorA       : Tcolor;                                                                    { bases }
    BaseColorT       : Tcolor;                                                                    { bases }
    BaseColorN       : Tcolor;                                                                    { bases }
    BaseColorGap     : Tcolor;                                                                    { bases }
    ColorSelectedB   : Tcolor;                                                                    { Selection }

    HighlightLowQV   : Boolean;                                                                   { Color bases with low QV in orange }
    HilightIUPACs    : Boolean;                                                                   { Highlight IUPAC bases in Pink }
    HilightMismat    : Boolean;                                                                   { Show mismatches in red in contig and in deep purple in reference. The input sequences will not be affected }

    {ASM RELATED}
    IsContig         : Boolean;                                                                   { True if this Cube is the output (contig) of the assembly process. False for all input sequences. }
    IsReference      : Boolean;                                                                   { True if this cube is a reference seq (need for: highlight reference in special color) }
    {VECINI}
    VeciniList       : TList;
    VeciniDepleated  : Boolean;                                                                   { arata daca am clasat toti vecinii acestui obiect }    { used by TAsmJobAbst.ClassifyObjectsRecursiv }
    ContigClass      : TBase;                                                                     { arata din ce contig face parte aceasta secventa }
    Assembled        : Boolean;                                                                   { True if this sequence was incorpored/assembled into the contig }
    {INIT}
    constructor Create(aLog: TRamLog);
    destructor  Destroy; override;
    procedure   Clear;   override;
    procedure   ClearSampleForReasembly;

    function GetDirection    (CONST DirectionChrPos: Integer): Char;                              { Determines the orientation of the sequence (Forward/Rev) based on its file name }

    {VECINI}
    function  vGetAllNeighbors: string;                                                           { intoarce toti vecinii acestui obiect, ca sir, separati cu virgula }
    function  vGetShortNameOf(VecinPoz: integer): string;
    procedure vAddVecin      (Vecin: TCubeObj);
    procedure vDeleteVecin   (PozVecinuluiInLista: integer);
    function  vGetPtr2Vecin  (IndexInLista: integer): TCubeObj;                                   { pointer catre un cub }

    {COLORS}
    procedure ComputeColors;                                                                      { Compute all colors for this object }
    procedure ComputeColorB(CONST BasePos: Integer);
    procedure ComputeColorT(CONST BasePos: Integer);
    property  ColorBkg     : TColor read FColorBkg write setColorBkg;                             { Color for normal cells (cells that are neither vectors (blue), erros (red), etc. This will give the overal color of the grid }
    property  ColorSolvedBkg: TColor read FSolvedBkg write setSolvedBkg;                          { Background color for 'solved' mismatches }

    function  CellHasLowQV(CONST BasePos: integer): Boolean;                                      { Returns true if the base has the QV under the 'EngTrim1.GoodQVTresh' threshold }
    property  RainbowTextClr : TColor  read FRainbowText   Write FRainbowText  default clOrange;  { Color of the text in rainbow mode }
    property  RainbowBkg     : Boolean read FRainbowBkg    Write FRainbowBkg   default FALSE;     { Cell background is colorful (it represents the ACGT colors) or it represents useful information (errors, mismatch points, low QV bases, etc) }
    property  ColorBookmarkT : TColor  read FBookmarkT     Write FBookmarkT    default clWhite;

    function  BuildQVGraph(ChromaWidth, ChromaHeight: Integer): TBitmap;                          { QV graph. Shows the quality of the sequence }

    {CONTIG RELATED}
    function  ContigNrOfBases: Integer;                                                           { The contig may contain dots. This count the bases without dots }
    function  ContigCount    : Integer;                                                           { Returns the number of subcontigs (if the contig contains multiple subcontigs) }
    function  ContigStarts   : Integer;                                                           { The contig may contain dots. Returns the position of the first base that is not dot }
    function  ContigEnds     : Integer;                                                           { The contig may contain dots. Returns the position of the last base that is not dot }
    function  ContigBases(InsertAds: Boolean)  : BaseString;                                                        { Returns the bases from which the contig is made, without the dots itmay hace at its ends. The dots in the middle (in case of multicontig) are not removed! }
    function  MismatchesPercent: Extended;
    function  MismatchesTotal: Integer;                                                           { Returns number of mismatches. The mismatches in the gray area (low quality ends) ARE ALSO counted. }
    function  Mismatches     : Integer;                                                           { Returns number of mismatches. The mismatches in the gray area (low quality ends) ARE NOT counted. }

    {ASM RELATED}
    function  AsmEnd         : Integer;                                                           { GridEnd= Grid.AsmOffset + NrOfBases. It refers to position in AsmJob; so the numeber returned represents an AsmGrid value not a Cub value }
    property  AsmOffset      : Integer read FAsmOffset     Write setAsmOffset  default 0;         { AsmGridOffset = Offset of this seq in the assembly grid/display. Indexed in 0. The TContigGrid will have to add 1 to this value, because on the position 0, it has the header. The AsmGrid don't have to add 1 because it is also indexed in 0. }
 end;



IMPLEMENTATION

uses CubeBase, cmMath, cGraphUtil, VectorProcessing;





{===============================================================================
                            CONSTRUCTOR / DESTRUCTOR
===============================================================================}
constructor TCubeObj.Create(alog: TRamLog);
begin
 inherited Create(aLog);                                                                          { about inherited: https://forums.codegear.com/thread.jspa?threadID=10046 }   { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }

 VeciniList     := TList.Create;
 FRainbowText   := clOrange;                                                                      { Color of the text in rainbow mode }
 BaseColorC     := CubicDNA.BaseColorC;                                                           {BASE COLOR}
 BaseColorG     := CubicDNA.BaseColorG;
 BaseColorA     := CubicDNA.BaseColorA;
 BaseColorT     := CubicDNA.BaseColorT;
 BaseColorN     := CubicDNA.BaseColorN;
 BaseColorGap   := CubicDNA.BaseColorGap;
 ColorEErrorB   := CellColorErrorB;                                                               { Error Color Bkg }
 ColorSolvedBkg := CellColorSolvedB;                                                              { Error Color Bkg }
 ColorSelectedB := CellColorSelectedB;                                                            { Selected }
 ColorBookmarkB := CellColorBookmarkB;
 ColorBookmarkT := CellColorBookmarkT;
 ColorBkg       := CellColorBkg;                                                                  { or is better clBtnFace ? }    { Color for normal cells (cells that are neither vectors (blue), erros (red), etc. This will give the overal color of the grid }
 FontColor      := CellFontColor;                                                                 { Font color for any text (other than ACGTN-) }
 FRainbowBkg    := FALSE;
 IsReference    := FALSE;
 IsContig       := FALSE;

 {OTHER COLORS}
 HighlightLowQV := TRUE;                                                                          { Color bases with low QV in orange }
 HilightMismat  := TRUE;
 HilightIUPACs  := TRUE;
end;

{ Apelez Clear cand:
     * creez un cub abstract (TCubeBase, TCubeImport, TCubObj), atunci apelez Clear imediat dupa Create.
     * creez un cubpersistent (TCubObjEx), atunci apelez Clear cand dau LoadFromFile sau Assign (daca o sa implementez vreodata Assign). }

destructor TCubeObj.Destroy;
begin
 FreeAndNil(VeciniList);
 inherited;
end;


procedure TCubeObj.Clear;
begin
 inherited Clear;

 { FOR GRID }
 FAsmOffset      := 0;                                                                            { deplasamentul in Grid }
 IsReference     := FALSE;
 { VECINI }
 Assembled       := FALSE;
 VeciniList.Clear;
 VeciniDepleated := FALSE;                                                                        { arata daca am clasat toti vecinii acestui obiect }
 ContigClass     := '?';                                                                          { arata din ce contig face parte aceasta secventa }
end;


procedure TCubeObj.ClearSampleForReasembly;                                                      {Nu pot sa apelez pru si simple 'Clear' caci asta ar face secventa invalida. Numai unele campuri trebuiesc resetate }
begin
 AsmOffset:= 0;
 VeciniDepleated:= FALSE;
 VeciniList.Clear;
 ContigClass:= '?';
 Assembled:= FALSE;

 { REMOVE ALL GAPS }                                                                              { Cris a zis ca trebuie sa scoata si GAP-urile inainte sa bage secventele din nou in asamblare din cauza ca gap-ul fute algoritmul de asamblare }
 RemoveAllGaps;
end;




function TCubeObj.getFileName: string;  { The contig does not have a 'FileName' because it is not an input file (it is not loaded from disk) }
begin
 Result:= inherited;
 {del  Verification
 if IsContig then RAISE Exception.Create('Contig name cannot be accessed via FileName! Use AsmJob.ContigFullName instead.');    }
end;






function TCubeObj.GetDirection(CONST DirectionChrPos: Integer): Char;                             { Determines the orientation of the sequence (Forward/Rev) based on its file name }
begin
 if DirectionChrPos= 0
 then raise Exception.Create('Invalid orientation parameters!');

 if DirectionChrPos <= Length(ScreenName)
 then
  begin
    Result:= upcase(ScreenName[DirectionChrPos]);
    if (Result<> 'F') AND (Result<> 'R')
    then Result:= '?';
  end
 else Result:= '?';
end;




{===============================================================================
                                 - VECINI -
================================================================================}
procedure TCubeObj.vAddVecin(Vecin: TCubeObj);
Begin
 VeciniList.Add(Vecin);                                                                           { ADD VECIN }
 //  Assembled:= (VeciniList.Count <= 0);                                                         { RAHAT !!!! Mai trebuie asta? se pare ca nu }
end;

procedure TCubeObj.vDeleteVecin(PozVecinuluiInLista: integer);
Begin
 VeciniList.Delete(PozVecinuluiInLista);
 //  Assembled:= (VeciniList.Count <= 0);
end;

function TCubeObj.vGetPtr2Vecin(IndexInLista: integer): TCubeObj;
Begin
  Result:= VeciniList[IndexInLista];
end;

function TCubeObj.vGetShortNameOf(VecinPoz: integer): string;                                     { Screen name }
begin
 Result:= TCubeObj(VeciniList[VecinPoz]).ScreenName;
end;

function TCubeObj.vGetAllNeighbors: string;
VAR i: Integer;
begin
 Result:= '';
 for i:= 0 TO VeciniList.Count-1 DO
  Result:= Result+ TCubeObj(VeciniList[i]).ScreenName+ ', ';

 if Length(Result)> 1
 then Result:= system.COPY(Result, 1, Length(Result)-2);                                                 { Remove the last character (which is ', ') }
end;








{--------------------------------------------------------------------------------------------------
                                 COMPUTE COLOR
--------------------------------------------------------------------------------------------------}
procedure TCubeObj.computeColorB(CONST BasePos: Integer);  { Background color }
begin
 if RainbowBkg
 then CellColorB[BasePos]:= GetBaseColor(Base[BasePos])
 else

  { BOOKMARK }
  if CellBookmark[BasePos]                                                                        { PRIORITATE MAXIMA }
  then CellColorB[BasePos]:= ColorBookmarkB
  else

    { VECTOR }
    if CellHasVector[BasePos] AND VectorShow
    then
      if CellTrimmed[BasePos]
      then CellColorB[BasePos]:= MixColors(VectorColor, CellColorDisabled, 90)                    { mix between vector and gray areas }
      else CellColorB[BasePos]:= VectorColor
    else

      { MISMATCH }
      if HilightMismat AND (CellMism[BasePos]= msMisSolved)
      then
       if CellTrimmed[BasePos]
       then CellColorB[BasePos]:= FColorESolvedGray                                              { combinatie intre solved si zone gri }
       else CellColorB[BasePos]:= ColorSolvedBkg
      else

       { MISMATCH UNSOLVED }
       if HilightMismat AND (CellMism[BasePos]= msMisUnsolved)
       then
         if isReference                                                                          { error in Reference }
         then CellColorB[BasePos]:= ReferenceSnpMims                                             { REFERENCE MUTATION ! }
         else CellColorB[BasePos]:= ColorEErrorB
       else

         { GRAY ENDS }
         if CellTrimmed[BasePos]
         then
           if isReference
           then CellColorB[BasePos]:= ReferenceLowQV                                             { combinatie intre referinta si zone gri }
           else CellColorB[BasePos]:= CellColorDisabled
         else

           { MUTATION SNP }
           if HilightIUPACs
           AND CellIsIUPAC(BasePos)
           then CellColorB[BasePos]:= CellColorSNP
           else

             { MUTATION CONTIG INDEL }
             if HilightIUPACs
             AND CellIsMutIndel(BasePos)
             then CellColorB[BasePos]:= CellClrMutIndel
             else

               { MUTATION CONTIG DIFFERENCE }
               if HilightIUPACs AND CellIsMutDiff(BasePos)
               then CellColorB[BasePos]:= CellClrMutDiff
               else

                { LOW QV }
                if CellHasLowQV(BasePos) AND HighlightLowQV
                then CellColorB[BasePos]:= CellColorLowQv
                else

                  { REFERENCE }
                  if IsReference
                  then CellColorB[BasePos]:= ReferenceCol
                  else
                                                                                                    { CONTIG }
                    if IsContig                                                                     { CONTIG }
                    then CellColorB[BasePos]:= FAutoContigColor
                    else CellColorB[BasePos]:= ColorBkg;                                            { Fundal }    { Culoarea celulei = culoarea Grid-ului }
end;


{ FONT COLOR }
procedure TCubeObj.computeColorT(CONST BasePos: Integer);
begin
 if CellBookmark[BasePos]
 then CellColorT[BasePos]:= ColorBookmarkT
 else
   if RainbowBkg
   then CellColorT[BasePos]:= RainbowTextClr
   else
     if HilightIUPACs
     AND CellIsMutant(BasePos)
     then CellColorT[BasePos]:= clBlack                                                           { Print SNP in black so it will be more visible on the pink background }
     else CellColorT[BasePos]:= GetBaseColor(Base[BasePos]);                                      { Normal colors }
end;



procedure TCubeObj.ComputeColors;                                                                 { Compute all colors for this object. Used in: TContigGrid.AssignCub, TContigGrid.ReLoadCub  }
VAR cl: Integer;
begin
 for cl:= 1 to NoOfBases DO
  begin
   ComputeColorB(cl);                                                                             { COMPUTE BKG COLOR }
   ComputeColorT(cl);                                                                             { COMPUTE FONT COLOR }
  end;
end;


procedure TCubeObj.setColorBkg(const Value: TColor);
VAR
   R,G,B: byte;
begin
 FColorBkg:= Value;

 { Make contig bkg color a bit bluish than standard color }  { This makes sense only if this cube is a contig }
 SplitColor2RGB(FColorBkg, R,G,B);

 if R < 12
 then R:= 0
 else R:= R- 12;

 if G < 12
 then G:= 0
 else G:= G- 12;

 if B > 243
 then B:= 255
 else B:= B+ 18;

 FAutoContigColor:= RGB(R,G,B);                                                                  {Culoare liniei contigului. E obtinuta automat din culoarea backgroundului Grid-ului }
end;


procedure TCubeObj.setSolvedBkg(const Value: TColor);
begin
 FSolvedBkg := Value;
 FColorESolvedGray:= DarkenColor(FSolvedBkg, 80);                                                {set this al 80% of FSolvedBkg }
end;


function TCubeObj.GetBaseColor(Base: TBase): TColor;                                             {Return the corresponding color for the A, C, G, T and N bases }
begin
 case UpCase(Base) of
   'C': Result:= BaseColorC;
   'G': Result:= BaseColorG;
   'A': Result:= BaseColorA;
   'T': Result:= BaseColorT;
   'N': Result:= BaseColorN;
   Gap: Result:= BaseColorGap;
 else
    Result:= clGray;                                                                             {orice alta baza care nu e C,G,A,T,GAP - de exemplu bazele ambigue cum ar fi W,S,M... }
 end;
end;








{--------------------------------------------------------------------------------------------------
   OTHERS
--------------------------------------------------------------------------------------------------}
function TCubeObj.CellHasLowQV(CONST BasePos: integer): Boolean;                                  { Returns true if the base has the QV under the 'EngTrim1.GoodQVTresh' threshold }
begin
 Result:=
  ((Base[BasePos]= 'N') AND NOT CellTrimmed[BasePos])
  OR
  (QVExist
     AND (CellQV[BasePos] < EngTrim1.GoodQVTresh)
     AND (NOT CellTrimmed[BasePos]));                                                             { daca celula asta e deja in zona Gri atunci nu are ros sa o mai pun cu rosu }
end;


function TCubeObj.AsmEnd: Integer;                                                                { Returns the position of the last base in AsmJob. It refers to position in AsmJob; so the number returned represents an AsmGrid value not a Cub value }
begin
 Result:= AsmOffset+ NoOfBases- 1;                                                                { -1 because in cub  bases are indexed in 1 but in AsmGrid they are indexed in 0 }
end;


procedure TCubeObj.setAsmOffset(Value: Integer);                                                  { AsmGridOffset = Offset in AsmGrid. Indexed in 0. The TContigGrid will have to add 1 to this value, because on the position 0, it has the header. The AsmGrid i don't have to add 1 because it is also indexed in 0. }
begin
 FAsmOffset:= Value;
end;











{--------------------------------------------------------------------------------------------------
   CONTIG
--------------------------------------------------------------------------------------------------}
function TCubeObj.ContigCount: Integer;                                                           { Returns the number of subcontigs (if the contig contains multiple subcontigs) }
VAR
   sBases, AllBases: BaseString;
   i: Integer;
begin
 Result:= 0;
 sBases:= '';
 AllBases:= Bases;

 for i:= 1 to Length(AllBases) DO
   begin
    { Build sub-contig }
    if AllBases[i]<> noCntgBase
    then sBases:= sBases+ AllBases[i];

    { Save sub-contig }
    if  (sBases> '')
    AND ((AllBases[i]= noCntgBase) OR (i= Length(AllBases)))
    then
     begin
       sBases:= '';
       inc(Result);
     end;
   end;
end;


{ Returns the contig.
    Gaps                  : removed
    Recognition sequences : removed (based on project settings)
    Dots at the end       : removed
    Dots in the middle    : not removed (multi contig case)
    Ads                   : inserted as specified }
function TCubeObj.ContigBases(InsertAds: Boolean): BaseString;
VAR i: Integer;
begin
 Assert(IsContig);
 Result:= '';

 for i:= ContigStarts to ContigEnds
   DO Result:= Result+ Base[i];

 { We remove the vectors }
 if Vectors.Detector<> NIL
 then Result:= TVectorDetector(Vectors.Detector).CleanVectors(Result);

 Result:= RemoveGaps(Result);

 { Ads }
 if InsertAds
 then Result:= InsertAdsOverBases(Result);                  { Insert ads if in DEMO mode }
end;



function TCubeObj.ContigNrOfBases: Integer;                                                       { The contig may contain dots. This count the bases without dots }
VAR i: Integer;
begin
 Assert(IsContig);

 Result:= 0;
 for i:= 1 to NoOfBases DO
   if base[i]<> noCntgBase
   then inc(Result);
end;


function TCubeObj.ContigStarts: Integer;                                                          { The contig may contain dots. Returns the position of the first base that is not dot }
VAR i: Integer;
begin
 Assert(IsContig);
 Result:= 0;

 for i:= 1 to NoOfBases DO
  if base[i]<> noCntgBase then
    begin
      Result:= i;
      Break;
   end;

 Assert(Result> 0);
end;


function TCubeObj.ContigEnds: Integer;                                                            { Returns the position of the last base that is not noCntgBase (dot). Indexed in 1. }
VAR i: Integer;
begin
 Assert(IsContig);
 Result:= 0;

 for i:= NoOfBases downto 1 DO
  if Base[i]<> noCntgBase then
    begin
      Result:= i;
      Break;
   end;

 Assert(Result> 0);
end;


function TCubeObj.MismatchesPercent: Extended;                                                    { Percent between number of good bases and ambiguities }
begin
 Result:= ProcentRepresent(Mismatches, NoOfBases);
end;


function TCubeObj.MismatchesTotal: Integer;                                                       { Returns number of mismatches. The mismatches in the gray area (low quality ends) ARE ALSO counted. }
VAR i: Integer;
begin
 Result:= 0;
 for i:= 1 to NoOfBases DO
   if CellMism[i]= msMisUnsolved
   then inc(Result);
end;


function TCubeObj.Mismatches: Integer;                                                            { Returns number of unsolved mismatches. The mismatches in the gray area (low quality ends) ARE NOT counted. }
VAR BasePos: Integer;
begin
 Result:= 0;
 for BasePos:= 1 TO NoOfBases-1 DO
 if NOT CellTrimmed[BasePos]                                                                      { nu calculez pentru cell-urile 'disabled/gri' }
    AND (CellMism[BasePos]= msMisUnsolved)                                                        { ErrorBkg }
 then inc(Result);
end;





{--------------------------------------------------------------------------------------------------
   QV Graph
--------------------------------------------------------------------------------------------------}
function TCubeObj.BuildQVGraph(ChromaWidth, ChromaHeight: Integer): TBitmap;
VAR x: Integer;
    CurBase, ScaleFactorX, ScaleFactorY: real;
begin
 Assert(HasChroma, 'TCubeObjEx.BuildQVGraph - Object has no chroma');
 Result:= TBitmap.Create;

 ScaleFactorX:= NoOfBases / ChromaWidth;
 ScaleFactorY:= 100 / ChromaHeight;
 Result.Width := ChromaWidth;
 Result.Height:= ChromaHeight;

 { BACKGROUND }
 Result.Canvas.Brush.Color:= TColor($D0FFD0); { lime }   { ChromaPasiveBkg; }
 Result.Canvas.FillRect(Rect(0, 0, ChromaWidth, ChromaHeight));

 { DRAW BAD QV BACKGROUND }
 Result.Canvas.Brush.Color:= BadEndsColor;
 Result.Canvas.FillRect(Rect(0, 0, round(GoodQVStart / ScaleFactorX), ChromaHeight));
 Result.Canvas.FillRect(Rect(round(LastGoodBase / ScaleFactorX), ChromaHeight, ChromaWidth, 0));

 { DRAW QV line }
 Result.Canvas.Pen.Color:= clNavy;
 Result.Canvas.MoveTo(0, ChromaHeight);
 x:= 0;
 CurBase:= 1;
 REPEAT
  Result.Canvas.LineTo(x, round(ChromaHeight- CellQV[round(CurBase)]/ScaleFactorY));
  inc(x);
  CurBase:= CurBase+ ScaleFactorX;
 UNTIL CurBase>= NoOfBases;
end;





END.(*==========================================================================



 TO IMPLEMENT IT:

   verific cazul Caroline in care ultimele 2 baze erau #0 si punctau catre  sample-ul 0 }

   iPointer:= SCF.BaseArray[BasePos].Ptr2Smpl;

   if iPointer< 1                                                                                 { verific cazul Caroline in care ulimile doau baze erau #0 si punctau catre  sample-ul 0 }
   then if  (BasePos-1> NoneAssigned)
        AND (SCF.BaseArray[BasePos-1].Ptr2Smpl+1< SCF.H.NrOfSamples)                              { point to the next sample after the last used sample }
        then iPointer:= SCF.BaseArray[BasePos-1].Ptr2Smpl+1
        else iPointer:= SCF.BaseArray[BasePos-1].Ptr2Smpl;
     -----------
   *)
