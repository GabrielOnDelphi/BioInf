UNIT ScfRead;

{===============================================================================
 Gabriel Moraru / Heracle BioSoft
 2016.08

 SCF READER
===============================================================================}

{todo 2: check comment generation in Fasta files }

INTERFACE

USES
   System.SysUtils, System.AnsiStrings, System.Classes, CubicDNA, ccStreamMem, scfbase;

TYPE
 TScfRead = class(TScfObj)
  private
  protected
    function  readHeader : Boolean;                               { Returns TRUE if magic number was found }
    procedure readSamples;                                        { Reads only SCF v3 samples }
    procedure readBases  ;                                        { Reads 'Bases' structures from a V3.x file }
    procedure deltaDeltaUnpack (VAR DiskTrace: TDiskTrace);       { Used by Reader }
  public
    function  LoadFromFile(CONST FullFileName: string): Boolean;
 end;




IMPLEMENTATION
USES ccBinary, ccCore, ccINIFile, clRamLog, FormLog;





{--------------------------------------------------------------------------------------------------
   Load From File
--------------------------------------------------------------------------------------------------}
function TScfRead.LoadFromFile(CONST FullFileName: string): Boolean;
begin
 Clear;
 Result:= FALSE;
 Assert(FileExists(FullFileName), 'SCF file does not exist: '+ FullFileName);

 { OPEN FILE }
 FFileName:= FullFileName;
 TRY
   FStream.LoadFromFile(FFileName);
   FStream.Position:=0;

   if FStream.Size< 512
   then Log.AddError(' Invalid file! The file is too small.')
   else
     { READ HEADER }
     if readHeader
     then
        if H.version < '3.00'
        then Log.AddError (' SCF version below 3.00 detected! Only SCF v3.00 or higher is supported.')
        else
          begin
            { READ SAMPLES }
            readSamples;

            { READ BASES }
            readBases;

            { READ Private Data }
            FStream.Position:= H.PrivateOffset;
            PrivateData:= string(FStream.ReadStringA(H.PrivateSize));

            { Done }
            Result:= NoOfBases > 2;
            if NOT Result
            then Log.AddError(' Sample file is empty!');
          end
     else
        { Maybe it is an ABI with wrong file ext? }
        if H.MagicNumber= 'abif'
        then Log.AddError(' This is an ABI file with wrong extension.' +CRLF+ 'Please correct your file name and try again.')
        else Log.AddError(' This is not a valid SCF file! Signature doesn''t match.');      { CONCLUSION: INVALID FILE }

 { Treat exception }
 if (NoOfBases<= 1)
 OR (NrOfSamples<= 1) then
  begin
    Log.AddError(' Invalid file: the sample is empty.');
    EXIT(FALSE);
  end;

 EXCEPT
   on E: Exception DO
    begin
     Log.AddError(' Invalid SCF file. Please send the file to us and we will try to recover it. '+ FullFileName);  { if I have an error, I put the file name at the beginning of the error message }
     Result:= FALSE;
    end;
 END;
end;



{--------------------------------------------------------------------------------------------------
   Read header
--------------------------------------------------------------------------------------------------}
function TScfRead.ReadHeader: Boolean;                                 { Return TRUE daca a gasit numarul magic}
CONST
   ScfMagic: AnsiString= '.scf';
VAR
   i: Integer;
begin
 FStream.Position:= 0;
 FStream.Read(H.MagicNumber, 4);
 Result:= SameText(H.MagicNumber, ScfMagic);

 if Result then
  begin                                                               // Some examples of real values:
   {# START OF CRITICAL READ ORDER }
   NrOfSamples       := FStream.RevReadInt;                           // 9674
   H.SamplesOffset   := FStream.RevReadInt;                           // 128
   NoOfBases         := FStream.RevReadInt;                                                      { BaseArray is Indexed in 1    -   This will make: SetLength(BaseArray, H.NoOfBases+ IndexedIn1) }
   H.ClipLeft        := FStream.RevReadInt;                           // 0                       { This info is not stored to cube }
   H.ClipRight       := FStream.RevReadInt;                           // 0                       { This info is not stored to cube }
   H.BasesOffset     := FStream.RevReadInt;                           // 77520
   H.CommentsSize    := FStream.RevReadInt;                           // 20
   H.CommentsOffset  := FStream.RevReadInt;                           // 86736
   FStream.Read (H.version, 4);                                       // 3.00
   H.SampleSize      := FStream.RevReadInt;                           // 2
   H.code_set        := FStream.RevReadInt;                           // 0
   H.PrivateSize     := FStream.RevReadInt;                           // 0
   H.PrivateOffset   := FStream.RevReadInt;                           // 0

   { Load spare data }
   for i:= 0 to 17
     DO H.Spare[i]:= FStream.RevReadInt;
   {# END OF CRITICAL READ ORDER }

   { Validate Sample Size }
   if (H.SampleSize < 1) OR (H.SampleSize >2) then
    begin
     Log.AddError(' SCF is malformed! Invalid samples size.');
     H.SampleSize:= 1;    { Fake things for older style SCF }
    end;

   { Read comments }
   if H.CommentsOffset+ H.CommentsSize <= Cardinal(FStream.Size)      { Make sure don't read outside the file }
   then                                                               { Had a case where the file was incomplete }
    begin
     FStream.Position:= H.CommentsOffset;
     Comments:= string(FStream.ReadStringA(H.CommentsSize));          { ReadStringU actually returns an ANSIString }
     Comments:= System.SysUtils.AdjustLineBreaks(Comments);
    end
   else
     Result:= FALSE;
 end;
end;



{--------------------------------------------------------------------------------------------------
   Read samples
--------------------------------------------------------------------------------------------------}
procedure TScfRead.readSamples;
Var i: Integer;
    DiskTraceA: TDiskTrace;                    { The 4 traces form together the Sample Matrix. Two bytes though in some cases 1 byte would be enough }
    DiskTraceC: TDiskTrace;
    DiskTraceG: TDiskTrace;
    DiskTraceT: TDiskTrace;
begin
 SetLength(DiskTraceA, NrOfSamples);           { Indexed in 0 }
 SetLength(DiskTraceC, NrOfSamples);
 SetLength(DiskTraceG, NrOfSamples);
 SetLength(DiskTraceT, NrOfSamples);

 {Position to where sample field starts }
 FStream.Position:= H.SamplesOffset;

 for i:= 0 TO NrOfSamples-1 DO
  begin
   FStream.Read( DiskTraceA[i], H.SampleSize); { Use ReadBuffer to read Count bytes from the stream into a buffer in cases where the number of bytes is known and fixed, for example when reading in structures. ReadBuffer is used internally for loading from a stream and copying from a stream. ReadBuffer calls Read to do the actual reading.  }
   DiskTraceA[i]:= Swap(DiskTraceA[i]);
  end;
 DeltaDeltaUnpack(DiskTraceA);

 for i:= 0 TO NrOfSamples-1 DO
  begin
   FStream.Read( DiskTraceC[i], H.SampleSize);
   DiskTracec[i]:= Swap(DiskTracec[i]);
  end;
 DeltaDeltaUnpack(DiskTraceC);

 for i:= 0 to NrOfSamples-1 DO
  begin
   FStream.Read( DiskTraceG[i], H.SampleSize);
   DiskTraceg[i]:= Swap(DiskTraceg[i]);
  end;
 DeltaDeltaUnpack(DiskTraceG);

 for i:= 0 to NrOfSamples-1 DO
  begin
   FStream.Read( DiskTraceT[i], H.SampleSize);
   DiskTracet[i]:= Swap(DiskTracet[i]);
  end;
 DeltaDeltaUnpack(DiskTraceT);

 { Transfer values }
 for i:= 0 to NrOfSamples-1 DO               { Transfer values from temp matrix (smallint) to final matrix (word) }
   begin
     TraceA[i]:= DiskTraceA[i];
     TraceC[i]:= DiskTracec[i];
     TraceG[i]:= DiskTraceg[i];
     TraceT[i]:= DiskTracet[i];
   end;
end;



{--------------------------------------------------------------------------------------------------
   Read Bases
   Read Bases structures from a V3.x file.
--------------------------------------------------------------------------------------------------}
procedure TScfRead.readBases;
VAR cr: AnsiChar;
    i, InvalidBases: Integer;
    LastPointer: Cardinal;
begin
 LastPointer := 0;
 InvalidBases:= 0;

 { Read pointers Bases->Samples }
 FStream.Position:= H.BasesOffset;
 for i:= IndexedIn1 to NoOfBases DO
  begin
   FStream.read ( BaseArray[i].Ptr2Smpl, 4);                                                  { RevReadInt }
   SwapCardinal ( BaseArray[i].Ptr2Smpl );                                                    { RevReadInt }

   { Treat exceptions: the case when the two bases point to the same sample }
   if  (LastPointer= BaseArray[i].Ptr2Smpl)
   AND (LastPointer < Cardinal(NrOfSamples) )                                                 { Make sure we don't exceed the number of samples}
   then BaseArray[i].Ptr2Smpl:= LastPointer+ 1;                                               { How we fix this: we let one of the bases to point to the original and the other one I make it point to the next sample }
   LastPointer:= BaseArray[i].Ptr2Smpl;
  end;

 { Read QVs from file }
 for i:= IndexedIn1 to NoOfBases DO FStream.read(BaseArray[i].prob_A, 1);
 for i:= IndexedIn1 to NoOfBases DO FStream.read(BaseArray[i].prob_C, 1);
 for i:= IndexedIn1 to NoOfBases DO FStream.read(BaseArray[i].prob_G, 1);
 for i:= IndexedIn1 to NoOfBases DO FStream.read(BaseArray[i].prob_T, 1);
 CheckQVs;

 { Read bases from file }
 for i:=IndexedIn1 to NoOfBases DO
   FStream.read(BaseArray[i].Base, 1);

 { Check the rest of the sequence for invalid bases }
 for i:= IndexedIn1 to NoOfBases DO
  begin
   cr:= BaseArray[i].Base;
   { Treat exceptions: We don't accept gaps in chromatogram. Replace them with 'N' }
   if NOT CharInSet(cr, DNA_Ambig) OR (cr= Gap) then
    begin
     BaseArray[i].Base:= 'N';
     inc(InvalidBases);
    end;
  end;

 if InvalidBases> 0
 then Log.AddWarn('   Invalid bases detected and automatically fixed.');

 { Spare data }
 for i:= IndexedIn1 to NoOfBases DO FStream.Read( BaseArray[i].spare[0], 1 );
 for i:= IndexedIn1 to NoOfBases DO FStream.Read( BaseArray[i].spare[1], 1 );
 for i:= IndexedIn1 to NoOfBases DO FStream.Read( BaseArray[i].spare[2], 1 );
end;





{--------------------------------------------------------------------------------------------------
   DELTA DELTA
--------------------------------------------------------------------------------------------------}
procedure TScfRead.DeltaDeltaUnpack(VAR DiskTrace: TDiskTrace);
var i: integer;
    PrevSample: Integer;
    Recover  : Integer;
begin
 PrevSample:= 0;
 for i:= 0 to NrOfSamples-1 do
   begin
     Recover := DiskTrace[i]+ PrevSample;
     { Treat exceptions: For some sequences (rarely) this value goues out of range }
     if (Recover>  32767) OR (Recover< -32768)
     then Recover:= 0;
     DiskTrace[i]:= Recover;
     PrevSample:= DiskTrace[i];
   end;
 PrevSample:= 0;
 for i:= 0 to NrOfSamples-1 do
   begin
     Recover := DiskTrace[i]+ PrevSample;
     if (Recover>  32767) OR (Recover< 0)   { In theory, from now on we should not have nagative values anymore }
     then Recover:= 0;
     DiskTrace[i]:= Recover;
     PrevSample:= DiskTrace[i];
   end;
end;


end.

