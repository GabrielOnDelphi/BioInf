
UNIT CubeStreamv4;

{==================================================================================================
 Heracle BioSoft SRL
 2016.10.31

 Object capabilities:
    + Read itself from stream (from files generated with Baser v2 to Baser v4)

==================================================================================================}

INTERFACE
USES
   System.SysUtils, System.Classes, System.Contnrs, Vcl.Graphics, Cube, ccStreamMem, clRamLog, ccRichLog;

TYPE
 RHeader= packed record
    MagicNumber     : string[16];                                                             { MagicNumber = 'DNA BaserV2 File' }
    Version         : string[16];
    SizeProps       : Cardinal;                                                               { Used in v4. Unused in v5 }
    SizeChroma      : Cardinal;
    SizeCellMx      : Cardinal;                                                               { It is the number of cells * size of a cell }
    Tag             : array[1..1024] of byte;
   end;

 TCubeObjEx4= class(TCubeObj)
  protected
    procedure ReadV4(aStream: TCubicMemStream; Hdr: RHeader);
 end;

CONST
   ctMagicNumber   = 'DNA Baser CubEx ';                                                      { 16 chars }
   ctVersion3      = '3.00.21         ';                                                      { 16 chars }
   ctVersion40     = '4.00.00         ';                                                      { 16 chars }
   ctVersion41     = '4.00.01         ';                                                      { 16 chars }    { Generated between Baser v4.1? and 4.36 (last v4 version) }
   ctCurrentVersion= '5.01.00         ';                                                      { 16 chars }    { This is the current version (v5.0). Bases are Unicode now }
  {ctDefChromaWidth= 50;    del
   ctDefRowHeight  = 16;    del }


IMPLEMENTATION
USES
   CubicDNA, CubeBase;

TYPE
 RCellV3 = record
    Base       : AnsiChar;                                                                         { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
    QV         : Byte;
    MismStatus : TMismStat;                                                                        { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
    Trimmed    : Boolean;
    OrigBase   : AnsiChar;
    EditPtr    : Integer;                                                                          { Stores the number of the edit that the user made in contig, edit which coresponds to this base (in imput sample) }
    Bookmark   : Boolean;
    Vector     : Boolean;                                                                          { True if this cell has a Vector assigned to it }
    ColorB     : TColor;                                                                           { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    ColorT     : TColor;                                                                           { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    Ptr2Smpl   : Integer;                                                                          { position of the sample associated with this base }
  end;

 RCell4 = record
    Base        : AnsiChar;                                                                        { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
    QV          : Byte;
    MismStatus  : TMismStat;                                                                       { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
    Trimmed     : Boolean;
    EdtOrig     : AnsiChar;                                                                        { Base before user edit. Makes sense for input samples only, not also for contig }
    EditPtr     : Integer;                                                                         { Makes sense for AsmJob. Points to the EditMX entry coresponding to the base edit (if this base was indeed edited) }
    Bookmark    : Boolean;
    Vector      : Boolean;                                                                         { True if this cell has a Vector assigned to it }
    ColorB      : TColor;                                                                          { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    ColorT      : TColor;                                                                          { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    Ptr2Smpl    : Integer;                                                                         { position of the sample associated with this base }
    SnpOrig_    : AnsiChar;                                                                        { If different than NoBase then it means that a SNP was detected here. This field does not store the SNP base but the original base as imported from SCF/ABI file. TConsensor.Conses will set this field to the same value as 'Base'  }
    CellProb    : RBaseProb;                                                                       { Probability of it being an A/C/G/T }
    SnpAreaRatio: Byte;                                                                            { Mutant peak to reference peak area ratio (%). I need it for the Log }
    SnpOverlap  : Byte;                                                                            { Overlap between mutant peak and reference peak (%). I need it for the Log  }
  end;

 RCell41 = record
    Base        : AnsiChar;                                                                        { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
    QV          : Byte;
    MismStatus  : TMismStat;                                                                       { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
    Trimmed     : Boolean;
    EdtOrig     : AnsiChar;                                                                        { Base before user edit. Makes sense for input samples only, not also for contig }
    EditPtr     : Integer;                                                                         { Makes sense for AsmJob. Points to the EditMX entry coresponding to the base edit (if this base was indeed edited) }
    Bookmark    : Boolean;
    Vector      : Boolean;                                                                         { True if this cell has a Vector assigned to it }
    ColorB      : TColor;                                                                          { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    ColorT      : TColor;                                                                          { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
    Ptr2Smpl    : Integer;                                                                         { position of the sample associated with this base }
    SnpOrig_    : AnsiChar;                                                                        { This field stores the original base as imported from SCF/ABI file. If different than NoBase then it means that a SNP was detected here.   DEL: TConsensor.Conses will set this field to the same value as 'Base'  }
    CtgMutation : TCtgMutation;  { NEW in v4.1 }
    CellProb    : RBaseProb;                                                                       { Probability of it being an A/C/G/T }
    SnpAreaRatio: Byte;                                                                            { Mutant peak to reference peak area ratio (%). I need it for the Log }
    SnpOverlap  : Byte;                                                                            { Overlap between mutant peak and reference peak (%). I need it for the Log  }
  end; 

 ACells3  = ARRAY of RCellV3;                                                                      { for DNA_Baser v3.0 cubes }
 ACells4  = ARRAY of RCell4;                                                                       { for DNA_Baser v4.0 cubes - Baser v4.0 to v4.8 }
 ACells41 = ARRAY of RCell41;                                                                      { for DNA_Baser v4.0 cubes - Baser v4.9 to v4.36 }


 RCubePropAnsi = packed record                                                                     { sizeof(RCubePropAnsi) is 978 }
   ParentType         : TBioFileType;
   ParentVersion      : string[32];
   OrigLength         : integer;
   Reversed           : Boolean;
   IsReference        : Boolean;
   IsPart             : Boolean;
   BasesNr            : Cardinal;
   EngTrim1           : REngTrim;
   EngTrim2           : REngTrim;
   NoOfSamples        : Cardinal;
   QVExist            : Boolean;
   GoodQVStart        : Integer;
   GoodQVEnd          : Integer;
   ContigClass        : AnsiChar;
   Assembled          : Boolean;
   Unused1            : Boolean;                                                                   { UNUSED. VeciniDepleated }
   AsmOffset          : Integer;
   IsContig           : Boolean;
   DateAdd            : TDateTime;                                                                 { UNUSED }
   DateModify         : TDateTime;                                                                 { UNUSED }
   OBSOLETE           : Boolean;
   Disabled           : Boolean;                                                                   { UNUSED }
   IsAssembled        : Boolean;
   Checked            : Boolean;                                                                   { UNUSED }
   xxxWrongString     : Cardinal;                                                                  { UNUSED. MetaData } { Here used to be a string. This is TOTALY WRONG hhere as I am not allowed to save/load from disk a dynamic string. Probably it reads 4 bytes }
   VectorsShow        : Boolean;                                                                   { show Vectors in blue }
   VectorsColor       : TColor;
   Legacy             : Boolean;                                                                   { In order to keep this record compatible with v3 }
   HighlightLowQV     : Boolean;
   RainbowBkg         : Boolean;
   HilightMismat      : Boolean;
   RainbowTextClr     : TColor;
   BaseColorC         : TColor;
   BaseColorG         : TColor;
   BaseColorA         : TColor;
   BaseColorT         : TColor;
   BaseColorN         : TColor;
   BaseColorGap       : TColor;
   ColorEErrorB       : TColor;
   ColorSolvedBkg     : TColor;
   ColorSelectedB     : TColor;
   ColorBookmarkB     : TColor;
   ColorBookmarkT     : TColor;
   ColorBkg           : TColor;
   FontColor          : TColor;
   {PAD}
   bTag1              : Boolean;
   bTag2              : Boolean;
   bTag3              : Boolean;
   iTag1              : Integer;
   iTag2              : Integer;
   iTag3              : Integer;
   iTag4              : Integer;
   iTag5              : Integer;
   iTag6              : Integer;
   iTag7              : Integer;
   sTag1              : string[255];  { This is written as: Counter + chars }
   sTag2              : string[255];
   sTag3              : string[255];
  end;



procedure TCubeObjEx4.ReadV4(aStream: TCubicMemStream; Hdr: RHeader);
VAR
   CubProp : RCubePropAnsi;
   CellMX3 : ACells3;
   CellMX40: ACells4;
   CellMX41: ACells41;
   i       : Integer;
begin
 { READ PROPERTIES }
 if sizeOf(CubProp) <> Hdr.SizeProps
 then RAISE Exception.Create('Cube stream size does not match record size!');

 aStream.Read(CubProp, Hdr.SizeProps);

 { READ DYNAMIC DATA }
 FileName  := String(aStream.ReadStringA);
 Comment   := String(aStream.ReadStringA);                                                    { COMMENTS }
 aStream.ReadStringA;                                                                         { ComentsEx. Unused }
 aStream.ReadStringA;                                                                         { Ads. Unused }
 aStream.ReadStringA;                                                                         { pad. Unused }
 RamLog.ImportRawData(String(aStream.ReadStringA));                                           { Log }
 aStream.ReadStringA;                                                                         { Unused }
 aStream.ReadStringA;                                                                         { Unused }

 NoOfBases:= CubProp.BasesNr;                                                                 { MEMORY ALOC }

 { READ CellMX v3}                                                                            { Transfer from v3 to current }
 if Hdr.Version= ctVersion3
 then
  begin
   SetLength(CellMX3, NoOfBases+ ctCellsIndex);
   aStream.ReadBuffer(CellMX3[0], Hdr.SizeCellMx);                                            { Load chroma. Am pus [0] pentru ca asa trebuie accesata o matrice dinamica. Daca ar fi fost statica nu as fiavut nevoie de asta. }
   for i:= ctCellsIndex to NoOfBases DO
    begin
     CellsMX[i].Base         := TBase(CellMX3[i].Base);                                       { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
     CellsMX[i].QV           := CellMX3[i].QV;
     CellsMX[i].MismStatus   := CellMX3[i].MismStatus;                                        { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
     CellsMX[i].Trimmed      := CellMX3[i].Trimmed;
     CellsMX[i].EdtOrig      := TBase(CellMX3[i].OrigBase);
     CellsMX[i].EditPtr      := CellMX3[i].EditPtr;                                           { Stores the number of the edit that the user made in contig, edit which coresponds to this base (in imput sample) }
     CellsMX[i].Bookmark     := CellMX3[i].Bookmark;
     CellsMX[i].Vector       := CellMX3[i].Vector;                                            { True if this cell has a Vector assigned to it }
     CellsMX[i].ColorB       := CellMX3[i].ColorB;                                            { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
     CellsMX[i].ColorT       := CellMX3[i].ColorT;                                            { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
     CellsMX[i].Ptr2Smpl     := CellMX3[i].Ptr2Smpl;                                          { position of the sample associated with this base }
     CellsMX[i].SnpOrig_     := noBase;
     CellsMX[i].CellProb.A   := 0;                                                            { Probability of it being an A }
     CellsMX[i].CellProb.C   := 0;                                                            { Probability of it being an C }
     CellsMX[i].CellProb.G   := 0;                                                            { Probability of it being an G }
     CellsMX[i].CellProb.T   := 0;                                                            { Probability of it being an G }
     CellsMX[i].SnpAreaRatio := 0;
     CellsMX[i].SnpOverlap   := 0;
   end;
  end
 else

 { READ CellMX v4.0 }                                                                         { Transfer from v4 to current }
 if Hdr.Version= ctVersion40
 then
  begin
   SetLength(CellMX40, NoOfBases+ ctCellsIndex);
   aStream.ReadBuffer(CellMX40[0], Hdr.SizeCellMx);                                           { Load chroma. Am pus [0] pentru ca asa trebuie accesata o matrice dinamica. Daca ar fi fost statica nu as fiavut nevoie de asta. }
   for i:= ctCellsIndex to NoOfBases DO
    begin
     CellsMX[i].Base         := TBase(CellMX40[i].Base);                                      { In mod normal obtin bazele din 'Base'. Doar can userul editeaza o baza manual, trec valoarea originala in 'OrigBase' si valoarea editata in 'Base' }
     CellsMX[i].QV           := CellMX40[i].QV;
     CellsMX[i].MismStatus   := CellMX40[i].MismStatus;                                       { color of the background: Red if the cell contains a mismatch, Green if the mismatch was fixed }
     CellsMX[i].Trimmed      := CellMX40[i].Trimmed;
     CellsMX[i].EdtOrig      := TBase(CellMX40[i].EdtOrig);
     CellsMX[i].EditPtr      := CellMX40[i].EditPtr;                                          { Stores the number of the edit that the user made in contig, edit which coresponds to this base (in imput sample) }
     CellsMX[i].Bookmark     := CellMX40[i].Bookmark;
     CellsMX[i].Vector       := CellMX40[i].Vector;                                           { True if this cell has a Vector assigned to it }
     CellsMX[i].ColorB       := CellMX40[i].ColorB;                                           { cand e Empty sau cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
     CellsMX[i].ColorT       := CellMX40[i].ColorT;                                           { cand userul nu vrea sa puna o culoare pentru fiecare baza CGAT in parte (RAINBOW MODE) }
     CellsMX[i].Ptr2Smpl     := CellMX40[i].Ptr2Smpl;                                         { position of the sample associated with this base }
     CellsMX[i].SnpOrig_      := TBase(CellMX40[i].SnpOrig_);
     CellsMX[i].CellProb     := CellMX40[i].CellProb;                                         { Probability of it being an ACGT }
     CellsMX[i].SnpAreaRatio := CellMX40[i].SnpAreaRatio;
     CellsMX[i].SnpOverlap   := CellMX40[i].SnpOverlap;

     if CellMX40[i].SnpOrig_<> noBase
     then CellsMX[i].CtgMutation:= cmDifference;
    end
  end
 else

 { READ CellMX v4.1-4.36 }                                                                   { Transfer from v4 to current }
 if Hdr.Version= ctVersion41
 then
  begin
   SetLength(CellMX41, NoOfBases+ ctCellsIndex);
   aStream.Read(CellMX41[0], Hdr.SizeCellMx);                                                 { Load chroma. Am pus [0] pentru ca asa trebuie accesata o matrice dinamica. Daca ar fi fost statica nu as fiavut nevoie de asta. }
   for i:= ctCellsIndex to NoOfBases DO
    begin
     CellsMX[i].Base         := TBase(CellMX41[i].Base);
     CellsMX[i].QV           := CellMX41[i].QV;
     CellsMX[i].MismStatus   := CellMX41[i].MismStatus;
     CellsMX[i].Trimmed      := CellMX41[i].Trimmed;
     CellsMX[i].EdtOrig      := TBase(CellMX41[i].EdtOrig);
     CellsMX[i].EditPtr      := CellMX41[i].EditPtr;
     CellsMX[i].Bookmark     := CellMX41[i].Bookmark;
     CellsMX[i].Vector       := CellMX41[i].Vector;
     CellsMX[i].ColorB       := CellMX41[i].ColorB;
     CellsMX[i].ColorT       := CellMX41[i].ColorT;
     CellsMX[i].Ptr2Smpl     := CellMX41[i].Ptr2Smpl;
     CellsMX[i].SnpOrig_      := TBase(CellMX41[i].SnpOrig_);
     CellsMX[i].CtgMutation  := CellMX41[i].CtgMutation;
     CellsMX[i].CellProb     := CellMX41[i].CellProb;
     CellsMX[i].SnpAreaRatio := CellMX41[i].SnpAreaRatio;
     CellsMX[i].SnpOverlap   := CellMX41[i].SnpOverlap;
    end
  end
 else
   RAISE Exception.Create('Unknown file version!');;

 { READ CHROMATOGRAM}
 NoOfSamples:= CubProp.NoOfSamples;                                                           { MEMORY ALOC }
 if Hdr.SizeChroma> 0
 then aStream.Read(Chroma[0], Hdr.SizeChroma);                                                { Load chroma. Am pus [0] pentru ca asa trebuie accesata o matrice dinamica. Daca ar fi fost statica nu as fiavut nevoie de asta. }                                  // Seek(F, Hdr.PrivateOffset);

 {= OBJECT PROPERTIES =}

 { SAMPLE }
 OrigLength     := CubProp.OrigLength;                                                        { Original length }
 ParentType     := CubProp.ParentType;
 Reversed       := CubProp.Reversed;                                                          { daca secventa a fost inversata }
 IsReference    := CubProp.isReference;
 IsPart         := CubProp.IsPart;
 { VECTORS }
 VectorShow     := CubProp.VectorsShow;                                                       { show Vectors in blue }
 VectorColor    := CubProp.VectorsColor;
 Vectors.Detector:= NIL;                                                                      { This will be restored later by assigning it to AsmJob.VectorDetector }
 { QV }
 EngTrim1       := CubProp.EngTrim1;
 EngTrim2       := CubProp.EngTrim2;
 FQVExist       := CubProp.QVExist;
 GoodQVStart    := CubProp.GoodQVStart;                                                       { This will make 'DirtyGoodBases:= TRUE' }
 LastGoodBase   := CubProp.GoodQVEnd;
 { BASES }
 DirtyBases     := TRUE;
 buildBases;                                                                                  { Rebuild bases }
 { CONTIG }
 ContigClass    := TBase(CubProp.ContigClass);                                                { arata din ce contig face parte aceasta secventa }
 Assembled      := CubProp.Assembled;
 AsmOffset      := CubProp.AsmOffset;                                                         { deplasamentul in Grid }
 IsContig       := CubProp.IsContig;
 { COLORS }
 HighlightLowQV := CubProp.HighlightLowQV;
 RainbowBkg     := CubProp.RainbowBkg;
 HilightMismat  := CubProp.HilightMismat;
 RainbowTextClr := CubProp.RainbowTextClr;
 BaseColorC     := CubProp.BaseColorC;
 BaseColorG     := CubProp.BaseColorG;
 BaseColorA     := CubProp.BaseColorA;
 BaseColorT     := CubProp.BaseColorT;
 BaseColorN     := CubProp.BaseColorN;
 BaseColorGap   := CubProp.BaseColorGap;
 ColorEErrorB   := CubProp.ColorEErrorB;
 ColorSolvedBkg := CubProp.ColorSolvedBkg;
 ColorSelectedB := CubProp.ColorSelectedB;
 ColorBookmarkB := CubProp.ColorBookmarkB;
 ColorBookmarkT := CubProp.ColorBookmarkT;
 ColorBkg       := CubProp.ColorBkg;
 FontColor      := CubProp.FontColor;
end;



end.
