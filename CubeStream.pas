
UNIT CubeStream;

{==================================================================================================
 Heracle BioSoft SRL
 2016.11.03

 Object capabilities:
    + Read/Write itself to stream  (from files generated with Baser v5)

 Object tree:
    * CubeBase -> CubeBaseQV -> CubeBaseSNP -> CubeImporter -> Cube -> CubeStream

==================================================================================================}

INTERFACE
USES
   Winapi.Windows, System.SysUtils, System.Classes, Cube, ccStreamMem, CubeStreamv4;

TYPE
 TCubeStream= class(TCubeObjEx4)
  private
    procedure ReadV5(aStream: TCubicMemStream; Hdr: RHeader);
  protected
  public
    procedure Save;
    procedure SaveToFolder (CONST DestFolder: string);
    procedure SaveToFile    (CONST FileName  : string);
    procedure WriteToStream (aStream: TCubicMemStream);
    function  ReadFromStream(aStream: TCubicMemStream): Boolean;
    function  LoadFromFile (CONST sFileName  : string): Boolean;
 end;


IMPLEMENTATION
USES
   CubicDNA, CubeBase, ccIO;




{--------------------------------------------------------------------------------------------------
   SAVE
--------------------------------------------------------------------------------------------------}
procedure TCubeStream.SaveToFile(CONST FileName: string);
VAR FStream: TCubicMemStream;
begin
 FStream:= TCubicMemStream.Create;
 TRY
   WriteToStream(FStream);
   FStream.SaveToFile(FileName);
 FINALLY
   FreeAndNil(FStream);
 END;
end;


   
procedure TCubeStream.Save;
begin
  SaveToFile(FileName);
end;



procedure TCubeStream.SaveToFolder(CONST DestFolder: string);
begin
 SaveToFile( Trail(DestFolder)+ ShortName );
end;



procedure TCubeStream.WriteToStream(aStream: TCubicMemStream);
VAR
   Hdr: RHeader;
begin
 Assert(Sizeof(rSample)= 12, 'CubeObjEx.WriteToStream - Invalid sample size!');
 Assert(aStream <> NIL, 'Stream is NIL');

 { HEADER }
 Hdr.MagicNumber  := ctMagicNumber;                                                { 16 chars }
 Hdr.Version      := ctCurrentVersion;                                             { 16 chars }
 Hdr.SizeChroma   := Length(Chroma) * Sizeof(rSample);                             { Documentatie: SizeOf knows if the data has been packed or not and returns the right value: http://www.delphibasics.co.uk/RTL.asp?Name=Packed }
 Hdr.SizeCellMx   := Length(CellsMX)* Sizeof(RCell5);
 Hdr.SizeProps    := 0;                                                            { Unused in v5 }
 FillMemory(@Hdr.Tag, SizeOf(Hdr.Tag), 0);                                         { Fill stuff with 0 so the record looks and compress better }
 aStream.Write(Hdr, SizeOf(Hdr));

 { FILE DATA}
 if iscontig
 then aStream.WriteStringU('Contig')     //del ExtractFileName(FFileName)       { The contig does not have a path and a file name. It depends on project's name/path }
 else aStream.WriteStringU(FileName);
 aStream.WriteStringU(Comment);                                                    { COMMENTS (from original ABI/SCF) }
 aStream.WriteInteger(OrigLength);                                                 { Original seq length, before trimming }
 aStream.WriteBoolean(Reversed);                                                   { True if the sequence was complement-reversed }
 aStream.WriteBoolean(IsReference);                                                { True if this cube is a reference seq (need for: highlight reference in special color) }
 aStream.WriteBoolean(IsPart);                                                     { True if this cube is a FASTA which is part of a multi-FASTA object }
 aStream.WriteByte(Ord(ParentType));                                               { What is the original provenience (input) of this Cube. Example: SCF, ABI, SEQ, FASTA, etc }

 { BASES }
 aStream.WriteInteger(NoOfBases);
 aStream.WriteInteger(EngTrim1.GoodBasesNr);                                       { Trim engine parameters }
 aStream.WriteInteger(EngTrim1.WindLength);                                        { Trim engine parameters }
 aStream.WriteInteger(EngTrim1.GoodQVTresh);                                       { Trim engine parameters }
 aStream.WriteInteger(EngTrim2.GoodBasesNr);                                       { Trim engine parameters }
 aStream.WriteInteger(EngTrim2.WindLength);                                        { Trim engine parameters }
 aStream.WriteInteger(EngTrim2.GoodQVTresh);                                       { Trim engine parameters }
 aStream.WriteBoolean(QVExist);                                                    { True if the original chromatogram has QV data }
 aStream.WriteInteger(GoodQVStart);                                                { The place where the first good base is located (gray/trimmed areas) }
 aStream.WriteInteger(LastGoodBase);                                               { The place where the last  good base is located }
 aStream.Write (CellsMX[0], Hdr.SizeCellMx);                                       { WRITE BASES }
 aStream.WriteCheckPoint;

 { CHROMA}
 aStream.WriteInteger(NoOfSamples);                                                { Number of chromatogram samples }
 
 { CONTIG}
 aStream.WriteChar(ansichar(ContigClass));
 aStream.WriteBoolean(Assembled);                                                  { True if this sequence was incorpored/assembled into the contig }
 aStream.WriteInteger(AsmOffset);                                                  { AsmGridOffset = Offset of this seq in the assembly grid/display. Indexed in 0. The TContigGrid will have to add 1 to this value, because on the position 0, it has the header. The AsmGrid don't have to add 1 because it is also indexed in 0. }
 aStream.WriteBoolean(IsContig);                                                   { True if this Cube is the output (contig) of the assembly process. False for all input sequences. }

 { VECTORS }
 aStream.WriteBoolean(VectorShow);                                                 { Show Vectors (in blue color) or not }
 aStream.WriteInteger(VectorColor);

 { COLORS }                                                                        { User color prefferences }
 aStream.WriteBoolean(HighlightLowQV);
 aStream.WriteBoolean(RainbowBkg);
 aStream.WriteBoolean(HilightMismat);
 aStream.WriteInteger(RainbowTextClr);
 aStream.WriteInteger(BaseColorC);
 aStream.WriteInteger(BaseColorG);
 aStream.WriteInteger(BaseColorA);
 aStream.WriteInteger(BaseColorT);
 aStream.WriteInteger(BaseColorN);
 aStream.WriteInteger(BaseColorGap);
 aStream.WriteInteger(ColorEErrorB);
 aStream.WriteInteger(ColorSolvedBkg);
 aStream.WriteInteger(ColorSelectedB);
 aStream.WriteInteger(ColorBookmarkB);
 aStream.WriteInteger(ColorBookmarkT);
 aStream.WriteInteger(ColorBkg);
 aStream.WriteInteger(FontColor);
 aStream.WriteInteger(FBookmarkT);
 aStream.WriteBoolean(HilightIUPACs);

 { WRITE CHROMA }
 if Hdr.SizeChroma> 0                                                              { SizeChroma = 0 means that this cube has no chromatogram assigned to it (comes from a FASTA file). The contig ALWAYS have a null chromatogram }
 then aStream.Write (Chroma[0], Hdr.SizeChroma);                                   { Save chroma. This is a dynamic array so we use [0] }

 { Check point }
 aStream.WriteCheckPoint;

 aStream.WritePadding(1024);
end;






{--------------------------------------------------------------------------------------------------
   LOAD
--------------------------------------------------------------------------------------------------}

function TCubeStream.LoadFromFile(CONST sFileName: string): Boolean;                 { Load TCubEx object from file }
VAR FStream: TCubicMemStream;
begin
 FileName:= sFileName;

 FStream:= TCubicMemStream.Create;
 TRY
   FStream.LoadFromFile(FileName);
   Result:= ReadFromStream(FStream);
 FINALLY
   FreeAndNil(FStream);
 END;
end;




function TCubeStream.ReadFromStream(aStream: TCubicMemStream): Boolean;
VAR Hdr: RHeader;
begin
 Result:= TRUE;

 { CHECK FILE SIZE }
 if aStream.Size < SizeOf(RHeader) then
  begin
   RamLog.AddError('Invalid file size!');
   EXIT(FALSE);
  end;

 { READ HEADER }
 aStream.Position:= 0;
 aStream.Read(Hdr, SizeOf(Hdr));

 { MAGIC NUMBER }
 if Hdr.MagicNumber <> ctMagicNumber then                                            { 16 chars }
  begin
   RamLog.AddError('Invalid magic number: '+ string(Hdr.MagicNumber));
   EXIT(FALSE);
  end;

 { FILE VERSION }
 if  (Hdr.Version<> ctCurrentVersion)                                                { current version }
 AND (Hdr.Version<> ctVersion40)                                                     { previous version }
 AND (Hdr.Version<> ctVersion41)                                                     { previous version }
 AND (Hdr.Version<> ctVersion3) then                                                 { previous version }
  begin
   RamLog.AddError('Unsupported version '+ string(Hdr.Version));
   EXIT(FALSE);
  end;

 { WHICH VERSION }
 if Hdr.Version= ctCurrentVersion                                                    { current version }
 then ReadV5(aStream, Hdr)
 else ReadV4(aStream, Hdr);
end;





procedure TCubeStream.ReadV5(aStream: TCubicMemStream; Hdr: RHeader);
begin
 Assert(Sizeof(rSample)= 12, 'CubeObjEx.WriteToStream - Invalid sample size!');

 { FILE DATA}
 FileName                := aStream.ReadStringU();
 Comment                 := aStream.ReadStringu;                                     { COMMENTS }
 OrigLength              := aStream.ReadInteger();
 Reversed                := aStream.ReadBoolean();
 IsReference             := aStream.ReadBoolean();
 IsPart                  := aStream.ReadBoolean();
 ParentType              := TBioFileType(aStream.ReadByte);

 { BASES }
 NoOfBases               := aStream.ReadInteger;
 EngTrim1.GoodBasesNr    := aStream.ReadInteger;                                     { QV DATA}
 EngTrim1.WindLength     := aStream.ReadInteger;
 EngTrim1.GoodQVTresh    := aStream.ReadInteger;
 EngTrim2.GoodBasesNr    := aStream.ReadInteger;
 EngTrim2.WindLength     := aStream.ReadInteger;
 EngTrim2.GoodQVTresh    := aStream.ReadInteger;
 FQVExist                := aStream.ReadBoolean;
 GoodQVStart             := aStream.ReadInteger;
 LastGoodBase            := aStream.ReadInteger;
 aStream.Read(CellsMX[0], Hdr.SizeCellMx);                                           { Read CELLSMX }
 aStream.ReadCheckPoint;
 DirtyBases              := TRUE;
 buildBases;                                                                         { Rebuild bases }

 { CHROMA}
 NoOfSamples             := aStream.ReadInteger;

 ContigClass             := char(aStream.ReadChar);                                  { CONTIG}
 Assembled               := aStream.ReadBoolean;
 AsmOffset               := aStream.ReadInteger;
 IsContig                := aStream.ReadBoolean;
 VectorShow              := aStream.ReadBoolean;                                     { VECTORS }
 VectorColor             := aStream.ReadInteger;

 { COLORS }
 HighlightLowQV          := aStream.ReadBoolean;
 RainbowBkg              := aStream.ReadBoolean;
 HilightMismat           := aStream.ReadBoolean;
 RainbowTextClr          := aStream.ReadInteger;
 BaseColorC              := aStream.ReadInteger;
 BaseColorG              := aStream.ReadInteger;
 BaseColorA              := aStream.ReadInteger;
 BaseColorT              := aStream.ReadInteger;
 BaseColorN              := aStream.ReadInteger;
 BaseColorGap            := aStream.ReadInteger;
 ColorEErrorB            := aStream.ReadInteger;
 ColorSolvedBkg          := aStream.ReadInteger;
 ColorSelectedB          := aStream.ReadInteger;
 ColorBookmarkB          := aStream.ReadInteger;
 ColorBookmarkT          := aStream.ReadInteger;
 ColorBkg                := aStream.ReadInteger;
 FontColor               := aStream.ReadInteger;
 FBookmarkT              := aStream.ReadInteger;
 HilightIUPACs           := aStream.ReadBoolean;  { Highlight IUPAC bases in Pink }


 { Read CHROMA }
 if Hdr.SizeChroma> 0                                                                { Contigs have a null chromatogram }
 then aStream.Read (Chroma[0], Hdr.SizeChroma);                                      { Save chroma. Am pus [0] pentru ca asa trebuie accesata o matrice dinamica. Daca ar fi fost statica nu as fiavut nevoie de asta.}

 { Check point }
 aStream.ReadCheckPoint;
 aStream.ReadPadding(1024);                                                          { Ad Unused }
end;




end.
