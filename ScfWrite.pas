
UNIT ScfWrite;

{===============================================================================
 Heracle BioSoft
 2016.07.13

 SCF WRITER
===============================================================================}

{DONE: Optimized SCF writer to work with Streams instead of direct file access }
       
INTERFACE

USES
   System.SysUtils, System.Classes,
   ccStreamMem, ccCore, ccINIFile, ScfBase;

TYPE
 TScfWriter = class(TScfObj)
  private
  protected
    procedure writeHeader;
    procedure writeSamples;
    procedure writeBases;
    procedure deltaDeltaCompact(VAR DiskTrace: TDiskTrace);                                        { used by Writer }
  public
    function SaveToFile(CONST FullFileName: string): Boolean;
 end;



IMPLEMENTATION
USES ccBinary;



{--------------------------------------------------------------------------------------------------
   CREATE
--------------------------------------------------------------------------------------------------}
function TScfWriter.SaveToFile(CONST FullFileName: string): Boolean;
CONST
   NonNulString= ' ';                                                                            { I can't write an empty string to disk so I have to insert a string }
begin
 FFileName:= FullFileName;

 if  (H.SampleSize<> 1)
 AND (H.SampleSize<> 2)
 then H.SampleSize:= 1;                                                                          { fake things for older style SCF }

 H.MagicNumber:= '.scf';
 Comments:= Comments+ NonNulString;                                                              { I can't write an empty string to disk so I have to insert a string }

 FStream.Position:= 0;
 writeHeader;
 writeSamples;
 writeBases;
 FStream.WriteChars(AnsiString(Comments));
 FStream.SaveToFile(FullFileName);

 Result:= TRUE;
end;



procedure TScfWriter.WriteHeader;
VAR DiskHdr: RHeaderSCF;
    i: integer;
begin
 Assert(SizeOf(H)= 128);

 { Copy the header }
 DiskHdr:= H;                                                                                      { because all bytes will be inversed (Little/Big Endian) I have to transfer the header in a new variable }
 DiskHdr.PrivateSize := 0;                                                                         { Nr. of bytes of Private data, 0 if none }
 DiskHdr.SamplesOffset  := SizeOf(DiskHdr);
 DiskHdr.BasesOffset    := DiskHdr.SamplesOffset  + H.NoOfSamples* DiskHdr.SampleSize* 4;
 DiskHdr.CommentsOffset := DiskHdr.BasesOffset    + H.NoOfBases  * SizeOf(TBase_v30);              { SCF Comments }
 DiskHdr.CommentsSize   := Length(Comments);
 DiskHdr.PrivateOffset  := DiskHdr.CommentsOffset+ DiskHdr.CommentsSize;                           { Byte offset from start of file}

 {Swap bytes}                                                                                      { From SCF version 3.0 the in memory structures and the data on the disk are not in the same format. See Little Endian vs Big Endian }
 {IN ACEST MOMENT INFORMATIA CONTINUTA IN DiskHdr NU MAI E VALIDA PT A FI FOLOSITA IN ACEASTA PROCEDURA !}
 SwapCardinal(DiskHdr.NoOfSamples);
 SwapCardinal(DiskHdr.SamplesOffset);
 SwapCardinal(DiskHdr.NoOfBases);
 SwapCardinal(DiskHdr.ClipLeft);                                                                   { This info was not stored to CUB. It is saved to disk even if it is meaningless }
 SwapCardinal(DiskHdr.ClipRight);                                                                  { This info was not stored to CUB. It is saved to disk even if it is meaningless } { Obsolete filed }
 SwapCardinal(DiskHdr.BasesOffset);
 SwapCardinal(DiskHdr.CommentsSize);
 SwapCardinal(DiskHdr.CommentsOffset);
 SwapCardinal(DiskHdr.SampleSize);
 SwapCardinal(DiskHdr.code_set);
 SwapCardinal(DiskHdr.PrivateSize);
 SwapCardinal(DiskHdr.PrivateOffset);

 for i:= 0 to 17
   DO SwapCardinal(DiskHdr.spare[i]);                                                              { I DON'T HAVE TO DO THAT: Unused }

 FStream.Write(DiskHdr, SizeOf(DiskHdr));
end;


procedure TScfWriter.WriteSamples;
Var i: Integer;
    DiskTraceA : TDiskTrace;                                                                       { 4 trace-rui formeaza un SAMPLE MATRIX } { pe 2 octeti desi in unele cazuri am nevoie doar de 1 byte }
    DiskTraceC : TDiskTrace;                                                                       { 4 trace-rui formeaza un SAMPLE MATRIX }
    DiskTraceG : TDiskTrace;                                                                       { 4 trace-rui formeaza un SAMPLE MATRIX }
    DiskTraceT : TDiskTrace;                                                                       { 4 trace-rui formeaza un SAMPLE MATRIX }

{sub}procedure CuanticWrite(Smpl: Smallint);                                                       { int2 = signed 16bit }
  begin
   if   H.SampleSize = 1
   then FStream.Write (Smpl, 1)                                                                    { Write  8 bit samples from a V3.x file }
   else
    begin
     Smpl:= System.Swap(Smpl);
     FStream.Write(Smpl, 2);                                                                       { Write 16 bit samples from a V3.x file }
    end;
  end;

begin
 SetLength(DiskTraceA, Length(TraceA));
 SetLength(DiskTraceC, Length(TraceC));
 SetLength(DiskTraceG, Length(TraceG));
 SetLength(DiskTraceT, Length(TraceT));

 { Transfer data from Trace to TraceADisk and convert from Word to SmallInt }
 for i:= 0 to NrOfSamples-1 DO                                                                     { Pastrez trace-urile originale nemodificate. Conversia DELTa se va face in TraceDisk }
   begin
    DiskTraceA[i]:= Smallint(TraceA[i]);
    DiskTraceC[i]:= Smallint(TraceC[i]);
    DiskTraceG[i]:= Smallint(TraceG[i]);
    DiskTraceT[i]:= Smallint(TraceT[i]);
   end;

 { TRACE A }
 DeltaDeltaCompact(DiskTraceA);
 for i:= 0 TO NrOfSamples-1
  DO CuanticWrite(DiskTraceA[i]);

 { TRACE C }
 DeltaDeltaCompact(DiskTraceC);
 for i:= 0 TO NrOfSamples-1
  DO CuanticWrite(DiskTraceC[i]);

 { TRACE G }
 DeltaDeltaCompact(DiskTraceG);
 for i:= 0 TO NrOfSamples-1
  DO CuanticWrite(DiskTraceG[i]);

 { TRACE T }
 DeltaDeltaCompact(DiskTraceT);
 for i:= 0 TO NrOfSamples-1
  DO CuanticWrite(DiskTraceT[i]);
end;


procedure TScfWriter.writeBases;
VAR i: integer;
begin
  for i:= IndexedIn1 to NoOfBases DO
   begin
    SwapCardinal(BaseArray[i].Ptr2Smpl);
    FStream.Write(BaseArray[i].Ptr2Smpl, 4);                                               { POINTERS }
   end;

  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].prob_A);                { QV A }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].prob_C);                { QV C }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].prob_G);                { QV G }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].prob_T);                { QV T }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteChar (BaseArray[i].Base  );                { BASE }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].spare[0]);              { PAD }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].spare[1]);              { PAD }
  for i:=IndexedIn1 to NoOfBases DO FStream.WriteByte (BaseArray[i].spare[2]);              { PAD }
end;





{--------------------------------------------------------------------------------------------------
   DELTA DELTA
--------------------------------------------------------------------------------------------------}
procedure TScfWriter.DeltaDeltaCompact(VAR DiskTrace: TDiskTrace);
var i: integer;
    LastSample: Integer;
    Encode  : Integer;
begin
 LastSample:= 0;
 for i:= 0 to NrOfSamples-1 do
   begin
     Encode := DiskTrace[i]- LastSample;
     if (Encode>  32767) OR (Encode< -32768)                                                       { Tratez cazul Heather. La unele secvente (foarte rar) imi iese din domeniu (Range error) asa ca trebuie sa fac verificarea asta }
     then Encode:= 0;
     LastSample:= DiskTrace[i];
     DiskTrace[i]:= Encode;
   end;

 LastSample:= 0;
 for i:= 0 to NrOfSamples-1 do
   begin
     Encode := DiskTrace[i]- LastSample;
     if (Encode>  32767) OR (Encode< -32768)                                                       { Tratez cazul Heather. La unele secvente (foarte rar) imi iese din domeniu (Range error) asa ca trebuie sa fac verificarea asta }
     then Encode:= 0;
     LastSample:= DiskTrace[i];
     DiskTrace[i]:= Encode;
   end;
end;




end.
