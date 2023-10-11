UNIT ReadABI; {new}

{=============================================================================================================
 Gabriel Moraru
 2016.07
==============================================================================================================


 Tester:
      c:\MyProjects\Biology\BASER\Resurse - Chromatogram tester\AbiTester.exe
      c:\MyProjects\Biology\BASER\Resurse - ABI format tester\AbiTester.exe
 Abi documentation:
      c:\MyProjects\Biology\Documentation\File formats\ABI\
 ==============================================================================}


INTERFACE
USES
   System.SysUtils, System.AnsiStrings, System.Classes, ccCore, ccINIFile, clRamLog, ccRichLog, ccStreamMem;

CONST
   ctCromaIndex = 1;                                                               { Chroma MX e indexata in 1 dar am de gand sa o indexez in zero }

TYPE
 { Directory }
 RDirEntry= packed record                                                          { Directory entry }{ 28 bytes. Old name= RAbiTag28 }
   Name      : string[4];                                                          { Obsolete } { Tag name }
   Number    : Longint;                                                            { Tag number. Signed 32 bits }
   ElemType  : Smallint;                                                           { Obsolete } { Element type code. Signed 16 bits }
   ElemSize  : Smallint;                                                           { Obsolete } { Size in bytes of an item }
   NrOfElems : Longint;                                                            { Number of directories. Original name: numelements }
   DataSize  : Longint;
   DataOffset: Longint;                                                            { Deplasamentul tag-ului in fisier. De obicei headerul se afla la sfarsitul fisierului }
   DataHandle: Longint;                                                            { Obsolete }
  end;


 { Header }
 TAbiHeader= record                                                                { The header has 128 bytes }
   MagicNumber    : string[4];                                                     { 4 cahracters ANSI string. Should contain "abif" }
   Version        : Word;
   MainDirectory  : RDirEntry;                                                     { 28 bytes }
   AbiTags        : string[94];                                                    { Array[1..94] of byte. I could not find documentation about it but most files seems to contain some tags in this space. in acest abi am gasit ca are 489 bytes: d:\Biology\SAMPLES\ Large sequence\ABI - Lopez - referinta de 10.000 baze\del\29 P103 17L c713.ab1 }
  end;


 { Trace }
 TAbiTrace= record                                                                 { There are 4 traces (AGCT) }
   Name         : Char;                                                            { WideChar is ok unless I save it to disk. Numele acestui trace: 'A', 'C', 'G', 'T' }
   NumberOf     : Longword;                                                        { Number of points - the bytes 8-11 represent a 32 bit integer containing the number of points in the trace value }
   Offset       : Longword;                                                        { Deplasamentul la care poate fi gasit acest camp }
   Points       : array of Word;                                                   { Array of sampling points. INDEXED IN 0 }
  end;


 { Side }
 TSide= record                                                                     { There are two sides: 1. Edited side &  2. Original side }
   NrOfBases      : Longword;                                                      { Number of bases (chars) in Bases matrix }  { unsigned 4-octeti integer}      { =Cardinal }
   NrOfPtrs       : Longword;                                                      { trebuie sa fie egal cu numarul de baze DECI asocierea se face de la BAZE CATRE SAMPLE-uri }
   BaseArrayOfset : Longword;                                                      { QvOffset pt matricea originala. Editarile sunt salvate in matricea 1}
   Ptr2Smpl       : array of Word;                                                 { Array of pointers from Bases to Samples (peaks) matrix. INDEXED in 1! }
   BaseArray      : AnsiString;                                                    { bazele }

   QvCount        : Longword;
   QvOffset       : Longword;
   QvMX           : array of Word;                                                 { Indexed in 0 }

   PtrsMXOffset   : Longword;                                                      { QvOffset for Ptr2Smpl matrix }
  end;



 TAbiObj = class(TObject)
  private
    FFileName: string;
    function  CorrectBaseOrder   (VAR Side: TSide): Boolean;
    function  CorrectPtr2SameSmpl(VAR Side: TSide): Boolean;                       { SEARCH FOR BASES POINTING TO THE SAME SAMPLE }
    procedure CheckPointersOutOfRange    (VAR Side: TSide);                                { VERIFIC SA NU DEPASEASCA NUMARUL SAMPLE-URILOR1 }
    procedure CheckFirstElem     (VAR Side: TSide);                                { VERIFIC PRIMUL ELEMENT}
    procedure CheckLastElem      (VAR Side: TSide);
    procedure readDirEntry       (VAR Entry: RDirEntry);
    procedure CheckPointerBaseEquality(Side: TSide);
  protected
    FStream: TCubicMemStream;
    function  showDirectoriesRawData: string;
    function  returnAddresses: string;                                             { For debugging }
    procedure writeBases   (Side: TSide; FileName: string);                        { For debugging }
  public
    H          : TAbiHeader;
    NoOfSamples: Longword;                                                         { Number of points - the bytes 8-11 represent a 32 bit integer containing the number of points in the trace value }
    Directories: array of RDirEntry;                                               { Aici se incarca TAG-urile propriu-zise }
    Trace1     : TAbiTrace;                                                        { Aici se incarca TRACE-urile propriuzise SI informatiile aferente }
    Trace2     : TAbiTrace;                                                        { 4 trace-rui formeaza un SAMPLE MATRIX }
    Trace3     : TAbiTrace;                                                        { 4 trace-rui formeaza un SAMPLE MATRIX }
    Trace4     : TAbiTrace;                                                        { 4 trace-rui formeaza un SAMPLE MATRIX }
    Side1      : TSide;                                                            { 1. Edited side }
    Side2      : TSide;                                                            { 2. Original side }
    Log        : TRamLog;
    constructor Create(aLog: TRamLog);                                             { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
    function  LoadFromFile(const FName: string; WriteLog: boolean): Boolean;       { Returns empty is file succesfully open or an error message if there were problems. Warnings are returned into the Warnings var. }
    procedure Clear;
    property  FileName: string read FFileName;
    function  AbiProperties: string;                                               { For debugging }
 end;


IMPLEMENTATION

USES
   System.Math, ccIO, CubicDNA, ccBinary;

CONST
   ctIrecuperabil= 'This ABI file is corrupt! Please send it to us and we will try to recover the information.';
   spatiere= 12;                                                                               { SPATIEREA medie intre doua baze. Asta ar trebuie calculata si nu pusa din burta! }




{===============================================================================
                            CONSTRUCTOR / DESTRUCTOR
===============================================================================}
constructor TAbiObj.Create;                                                                    { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
begin
 inherited Create;                                                                             { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }
 Log:= aLog;
 Clear;                                                                                        { asta e initializata ultima }
end;



procedure TAbiObj.Clear;

  procedure ClearSide(VAR Side: TSide);
  begin
   Side.NrOfBases:= 0;
   Side.NrOfPtrs := 0;
   Side.BaseArrayOfset:= 0;
   SetLength(Side.Ptr2Smpl, 0);
   Side.BaseArray:= '';
   Side.QvCount  := 0;
   Side.QvOffset := 0;
   SetLength(Side.QvMX, 0);
  end;

  procedure ClearTrace(VAR Trace: TAbiTrace);
  begin
    Trace.Name:= '?';
    Trace.NumberOf:= 0;
    Trace.Offset:= 0;
    SetLength(Trace.Points, 0);
  end;

begin
 FFileName   := '';
 NoOfSamples := 0;                                                                      { Number of points - the bytes 8-11 represent a 32 bit integer containing the number of points in the trace value }

 FillChar  (H, SizeOf(H), 0);
 SetLength (Directories, 0);

 ClearTrace (Trace1);
 ClearTrace (Trace2);
 ClearTrace (Trace3);
 ClearTrace (Trace4);

 ClearSide  (Side1);
 ClearSide  (Side2);
end;






{===============================================================================
                                    LOAD
===============================================================================}

function TAbiObj.LoadFromFile(const FName: string; WriteLog: boolean): boolean;                    { Returns empty is file succesfully open or an error message if there were problems. Warnings are returned into the Warnings var. }
VAR
    i: Integer;
    iDataOffset: Integer;

 function ReadPtrs(VAR Side: TSide): Boolean;
 VAR p: Integer;
 begin
  Result:= TRUE;

  if Side.NrOfPtrs< 1 then
   begin
    Log.AddError('Number of pointers to samples is: zero!');
    EXIT(FALSE);
   end;

  FStream.Position:= Side.PtrsMXOffset+ 1;
  SetLength(Side.Ptr2Smpl, Side.NrOfPtrs+ ctCromaIndex);

  { Tratez cazul 'Boulahtouf Abdelhay': am grija sa nu citesc in afara fisierului }
  p:= Side.PtrsMXOffset+ (Side.NrOfPtrs* 2);
  if p >= FStream.Size
  then Log.AddVerb('PtrsMXOffset beyond end of file! The program will attempt to recover data from file.');
  While (p >= FStream.Size) DO
   begin
    Dec(Side.NrOfPtrs);                                                                       { Decrease number of pointers until I am back inside the file }
    p:= Side.PtrsMXOffset+ (Side.NrOfPtrs* 2);
   end;

  { READ POINTERS }
  for p:= 1 to Side.NrOfPtrs DO
    FStream.Read (Side.Ptr2Smpl[p], 2);                                                       { NORMAL READ }
 end;


 procedure ReadTraces(VAR Trace: TAbiTrace);
 VAR t: Integer;
 begin
  SetLength(Trace.Points, Trace.NumberOf);
  //Unnecessary:
  //       We MUST fille the arrays with zeros because in some cases not all traces have the same number of elements. For details, see how NoOfSamples is initialized (to the longest trace): http://stackoverflow.com/questions/39960356/how-to-create-a-procedure-like-setlength-that-also-zeros-the-memory/39960664#39960664
  //       FillChar(Trace.Points[0], length(Trace.Points)* SizeOf(word), 0);

  FStream.Position:= Trace.Offset;                                                              { READ TRACE 1 }
  for t:= 0 TO Trace.NumberOf-1 DO
   begin
    FStream.Read(Trace.Points[t], 2);
    Trace.Points[t]:= System.Swap(Trace.Points[t]);
   end;
 end;


 procedure ReadQvMatrix(VAR Side: TSide);
 VAR t: Integer;
     QVReadSize: Cardinal;
 begin
  Assert(Side.QvOffset < FStream.Size- Side.NrOfBases, 'QV matrix out of range!');

  FStream.Position:= Side.QvOffset;

  { Do I have QVs? }
  if Side.QvCount = 0 then                                                          { in cazul in care campul PCON lipseste atunci fortez/incarc matricea de QV-uri cu zerouri - CAZUL Michi }
    begin
     Side.QvCount:= Side.NrOfBases;
     Log.AddInfo('This file does not contain information about quality of the called bases. It is HIGHLY recommended to activate the internal base caller!!');
    end;

  { QvCount = NoOfBases? }
  if Side.NrOfBases <> Side.QvCount
  then Log.AddVerb('QV matrix <> number of bases: '+ IntToStr(Side.NrOfBases)+ '. QVs: '+ IntToStr(Side.QvCount));

  { The size of QvMX MUST always be the same as the size of BasesMX }
  SetLength(Side.QvMX, Side.NrOfBases);   { If I have more QV than bases, I ignore the extra QVs. If I don't have enough QVs, I simply read what I have and the rest of the matrix will remained as it is (filled with zeros) }

  { Read data }
  //I tried to do this but won't work: FStream.Read (Side.QvMX[0], QVReadSize* SizeOf(Word));
  QVReadSize:= min(Side.NrOfBases, Side.QvCount);                         { I cannot read more than Side.QvCount because there is not more data on disk but I also cannot read more than NoOfBases because the size of the QVMatrix is as big as NoOfBases }
  for t:= 0 to QVReadSize-1
   DO FStream.Read (Side.QvMX[t], 1);
 end;


begin
 Clear;
 Result:= FALSE;
 FFileName:= FName;
 Assert(FileExists(FFileName));

 FStream:= TCubicMemStream.Create;
 TRY
 TRY

  FStream.LoadFromFile(FFileName);
  FStream.Position:= 0;
  if FStream.Size< 512 then     { File is too small? }                                                 { Size: Indicates the size in bytes of the stream }
   begin
    Log.AddError('The size of this ABI file is invalid!');
    EXIT;
   end;


  { === MAGIC NO =================================== }
  H.MagicNumber:= FStream.ReadStringA(4);
  if LowerCase(H.MagicNumber) <> 'abif' then       { Check valid MagicNumber }
   begin
    if H.MagicNumber= '.scf'                                                                           { verifica daca e un SCF cu extensie gresita }
    then Log.AddError('This is a SCF file with a wrong extension.'+ CRLF+ 'Please correct file''s extension and try again.')
    else Log.AddError(ctIrecuperabil);
    EXIT;
   end;


  { === READ VERSION =================================== }
  H.Version:= FStream.RevReadWord;


  { === READ MainDirectory ============================= }
  ReadDirEntry(H.MainDirectory);

  { Tratez cazul Caroline MPI }
  if H.MainDirectory.DataOffset >= FStream.Size then                                                   { am pus 'Cardinal' aici ca sa nu mai imi dea acel warning: Comparing blablalbla }
   begin
    Log.AddError('Premature end of file. '+ ctIrecuperabil);
    Exit;
   end;

  { In some files this will contain some text data from the machine }
  FStream.Position:= 129;
  H.AbiTags:= FStream.ReadStringA(94);
  H.AbiTags:= AnsiString(RemoveLowChars(string(H.AbiTags)));


  { === READ DIRECTORIES =============================== }
  SetLength(Directories, H.MainDirectory.NrOfElems);
  FStream.Position:= H.MainDirectory.DataOffset;

  for i:= 0 TO H.MainDirectory.NrOfElems-1 DO
   begin
    ReadDirEntry(Directories[i]);
    iDataOffset:= Directories[i].DataOffset;

    { Decode: NAME ASSOCIATED WITH EACH TRACE }
    if Directories[i].Name= 'FWO_' then                                                                { LOOKING FOR THE TYPE/ORDER OF TRACES - FWO field }
     begin
      Trace1.Name:= Char(GetByte(1, iDataOffset));
      Trace2.Name:= Char(GetByte(2, iDataOffset));
      Trace3.Name:= Char(GetByte(3, iDataOffset));
      Trace4.Name:= Char(GetByte(4, iDataOffset));
     end;

    { Decode: BASES INFO }
    if Directories[i].Name= 'PBAS' then
     if Directories[i].Number=1 then
       begin
         Side1.BaseArrayOfset := iDataOffset;                                                           { StadenName = BaseEntry }
         Side1.NrOfBases      := Directories[i].NrOfElems;
       end else
       begin
         Side2.BaseArrayOfset := Directories[i].DataSize;
         Side2.NrOfBases      := Directories[i].NrOfElems;
       end;

    { Decode: POINTERS INFO }
    if Directories[i].Name= 'PLOC' then                                                                { The PLOC field appears twice }
      if Directories[i].Number=1 then
      begin
        Side1.PtrsMXOffset  := iDataOffset;                                                            { PLOC offsetul la care se afla matricea de pointeri catre SAMPLEs }
        Side1.NrOfPtrs := Directories[i].NrOfElems;
       end
      else
       begin
        Side2.PtrsMXOffset := iDataOffset;
        Side2.NrOfPtrs := Directories[i].NrOfElems;
      end;

    { Decode: QV (CONfidence) }
    if Directories[i].Name= 'PCON' then                                                                { am gasit un caz in care masina nu a generat acest camp - CAZUL Michi }
      if Directories[i].Number=1 then
       begin
        Side1.QvOffset:= iDataOffset;
        Side1.QvCount:= Directories[i].NrOfElems;
       end
      else
       begin
        Side2.QvOffset:= iDataOffset;
        Side2.QvCount:= Directories[i].NrOfElems;
       end;

    { Decode: TRACES INFO }
    if Directories[i].Name= 'DATA' then
      if Directories[i].Number=9  then                                                                 { DATA segments 9 - 12 contains the adresses of traces. We don't know which segments represents which base. To know that you also have to look for the 28 bytes segments starting with "FWO_". The address of each segment is given as 32 bit long integer presnt form bytes 20-23 of each 28 byte segment. Note all the five offset addresses. In each 28 bytes of DATA segment the bytes 8-11 reprenstn a 32 bit integer containing the number of point in the trace value. }
        begin
         Trace1.Offset  := iDataOffset;
         Trace1.NumberOf:= Directories[i].NrOfElems;
        end else                                                                                       { the address of sample traces }
      if Directories[i].Number=10 then
        begin
         Trace2.Offset  := iDataOffset;
         Trace2.NumberOf:= Directories[i].NrOfElems;
        end else                                                                                       { the address of sample traces }
      if Directories[i].Number=11 then
        begin
         Trace3.Offset  := iDataOffset;
         Trace3.NumberOf:= Directories[i].NrOfElems;
        end else                                                                                       { the address of sample traces }
      if Directories[i].Number=12 then
        begin
         Trace4.Offset  := iDataOffset;
         Trace4.NumberOf:= Directories[i].NrOfElems;
        end;                                                                                           { the address of sample traces }
  end;


  { === NO OF SAMPLES =================================== }
  NoOfSamples:= max(max(Trace1.NumberOf, Trace1.NumberOf), max(Trace3.NumberOf, Trace4.NumberOf));     { SET HIGHEST NoOfSamples }                                                                { Trace1 should be = Trace2 = Trace3 = Trace4... }

  if (NoOfSamples <> Trace1.NumberOf)
  OR (NoOfSamples <> Trace2.NumberOf)
  OR (NoOfSamples <> Trace3.NumberOf)
  OR (NoOfSamples <> Trace4.NumberOf)
  then Log.AddVerb('Mismatch between the number of samples in the 4 traces! Expanding all other 3 traces to the length of the longest trace.');

  if NoOfSamples <= 5 then                                                                              { CHECK SAMPLES }
   begin
    Log.AddError('This ABI file has no chromatogram!');
    EXIT;
   end;


  { === BASES === }

  { Check no of bases }
  Result:= (Side1.NrOfBases >= ctMinimimSequence) AND (Side2.NrOfBases >= ctMinimimSequence);
  if NOT Result
  then Log.AddError('This ABI file has less less than '+ i2s(ctMinimimSequence) +' bases!');

  { Read bases }
  FStream.Position:= Side1.BaseArrayOfset;
  Side1.BaseArray := FStream.ReadStringA(Side1.NrOfBases);

  { Read bases }
  FStream.Position:= Side2.BaseArrayOfset;
  Side2.BaseArray := FStream.ReadStringA(Side2.NrOfBases);


  { === PTRS 2 SAMPLE ==== }
  if NOT ReadPtrs(Side1) then EXIT;
  if NOT ReadPtrs(Side2) then EXIT;


  { CHECKS }
  CheckPointerBaseEquality(Side1);                                                     { No of pointers = No of bases? }
  CheckPointerBaseEquality(Side2);

  CheckPointersOutOfRange(Side1);                                                      { VERIFIC SA NU DEPASEASCA NUMARUL SAMPLE-URILOR }
  CheckPointersOutOfRange(Side2);

  if CorrectPtr2SameSmpl(Side1)                                                        { SEARCH FOR BASES POINTING TO THE SAME SAMPLE }      { Trateaza cazul Assem in care doua baze din ABI.BaseArray puncteaza catre acelasi sample. ABI are pointeri de la Bases catre Samples insa eu am de la Samples catre Bases, asa ca data viitoare cand scanez matricea de Samples ca sa gasesc pointerii catre baze, o sa sar peste una din baze, caci Sample-ule meu are doar un pointer si nu poate sa puncteze catre ambele baze in acelasi timp. }
  OR CorrectPtr2SameSmpl(Side2)
  then Log.AddVerb('Bases pointing to the same sample found. Fixed.');

  CheckFirstElem(Side1);                                                               {VERIFIC PRIMUL ELEMENT}
  CheckFirstElem(Side2);

  CorrectBaseOrder(Side1);                                                             { CORECTEZ ORDINEA BAZELOR }
  CorrectBaseOrder(Side2);

  CheckLastElem(Side1);                                                                {VERIFIC ULTIMUL ELEMENT}
  CheckLastElem(Side2);


  { === QV MATRIX === }
  if (Side1.QvOffset > 0) AND (Side1.QvCount > 0)
  then ReadQvMatrix(Side1);
  if (Side2.QvOffset > 0) AND (Side2.QvCount > 0)
  then ReadQvMatrix(Side2);


  { === TRACES === }
  ReadTraces(Trace1);
  ReadTraces(Trace2);
  ReadTraces(Trace3);
  ReadTraces(Trace4);


  Log.AddHint('We strongly recommend you to use the SCF format instead of ABI!');
 FINALLY
   FreeAndNil(FStream);
 END;
 EXCEPT
   on E: Exception DO
    begin
     Log.AddError('Invalid ABI file. Please send this file to us and we will try to recover the information from it.');                                { Excpetion is handled here and does not propagate anymore }
     Result:= FALSE;
    end;
 END;
end;



procedure TAbiObj.ReadDirEntry(VAR Entry: RDirEntry);
begin
  Entry.Name      := ShortString(FStream.ReadStringA(4));
  Entry.Number    := FStream.RevReadLongInt;
  Entry.ElemType  := FStream.RevReadSmallInt;
  Entry.ElemSize  := FStream.RevReadSmallInt;
  Entry.NrOfElems := FStream.RevReadLongInt;
  Entry.DataSize  := FStream.RevReadLongInt;
  Entry.DataOffset:= FStream.RevReadLongInt;
  Entry.DataHandle:= FStream.RevReadLongInt;
//  EnsureByte()
end;









{===============================================================================
   VERIFICATIONS
===============================================================================}

{ 1. }
procedure TAbiObj.CheckPointerBaseEquality(Side: TSide);                    { No of pointers = No of bases? }
begin
  if Side.NrOfBases < side.NrOfPtrs
  then
   begin
    Log.AddVerb('No of bases ('+ IntToStr(Side.NrOfBases)+') <> No of pointers  ('+ IntToStr(Side.NrOfPtrs)+') !');
    Side.NrOfPtrs:= Side.NrOfBases;
   end
  else
    if Side.NrOfBases > side.NrOfPtrs
    then Log.AddVerb('!No of bases ('+ IntToStr(Side.NrOfBases)+') <> No of pointers  ('+ IntToStr(Side.NrOfPtrs)+') ! (Side)');
end;


{ 2. }
procedure TAbiObj.CheckPointersOutOfRange(VAR Side: TSide);                                               { Fix pointers pointing to invalid samples }
VAR median: Integer;
    p: Cardinal;
begin
 { FIRST we need to check for negative pointers }
 for p:= 1 to Side.NrOfPtrs DO   {NOTE: original code had here   NoOfSamples and it worked! }
  if Side.Ptr2Smpl[p]< ctCromaIndex then
   begin
    Log.AddVerb(Format('Warning in %s file! The base %s on position %d points to sample 0', [ExtractFileName(FileName), Side.BaseArray[p], p]));
    Side.Ptr2Smpl[p]:= ctCromaIndex;
   end;

 { Some ABIs points to a sample that does not exists (higher than NoOfSamples) }
 for p:= 1 to Side.NrOfPtrs DO
  if Side.Ptr2Smpl[p]> NoOfSamples                                  { jeffrey reported a sample 'p61_1_CP43_SP1F_O14.ab1' that has this problem }
  then
    begin
      Log.AddVerb(Format('The base %d points to sample %d which is out of range. Current number of samples %d. Fixed.', [p, Side.Ptr2Smpl[p], NoOfSamples]));

      if  (p> 1)                                   { Because I use p-1 below }
      AND (p< Side.NrOfPtrs)                       { Because I use p+1 below }
      then median:= abs(Side.Ptr2Smpl[p+1] - Side.Ptr2Smpl[p-1])
      else median:= spatiere;

      Side.Ptr2Smpl[p]:= Side.Ptr2Smpl[p-1]+ (median DIV 2);
    end;
end;



{ 3. SEARCH FOR BASES POINTING TO THE SAME SAMPLE }
function TAbiObj.CorrectPtr2SameSmpl(VAR Side: TSide): Boolean;      { Trateaza cazul Assem in care doua baze din ABI.BaseArray puncteaza catre acelasi sample. ABI are pointeri de la Bases catre Samples insa eu am de la Samples catre Bases, asa ca data viitoare cand scanez matricea de Samples ca sa gasesc pointerii catre baze, o sa sar peste una din baze, caci Sample-ule meu are doar un pointer si nu poate sa puncteze catre ambele baze in acelasi timp. }
VAR
   CurrentPointer: Word;
   NoErrors: Boolean;
   nextBase, curBase: Cardinal;
begin
 Result:= FALSE;
 REPEAT
  NoErrors:= TRUE;
  for CurBase:= 1 to Side.NrOfBases-1 DO
   begin
    NextBase:= CurBase+ 1;
    CurrentPointer:= Side.Ptr2Smpl[CurBase];

    { Scan all remaining BASES to see if another one points to the same SAMPLE }
    WHILE (NextBase<= Side.NrOfBases)
      AND (Side.Ptr2Smpl[NextBase]<> 0)
      AND (Side.Ptr2Smpl[NextBase]= CurrentPointer) DO
      begin
        Result:= TRUE;
        NoErrors:= FALSE;
        if CurrentPointer+1<= NoOfSamples                                                          { make sure you don't assign to a sample that doesn't exists }
        then Side.Ptr2Smpl[NextBase]:= CurrentPointer+ 1                                           { make this BASE to point to the next SAMPLE }
        else Side.Ptr2Smpl[NextBase]:= 0;
        Log.AddVerb('Invalid pointer to sample: '+ IntToStr(Side.Ptr2Smpl[NextBase]) +'. Fixed.');
        Inc(NextBase);
      end;
   end;
 UNTIL NoErrors;
end;




procedure TAbiObj.CheckFirstElem(VAR Side: TSide);                                               {VERIFIC PRIMUL ELEMENT}
VAR i: Integer;
begin
 if (Side.Ptr2Smpl[ctCromaIndex+1]- Side.Ptr2Smpl[ctCromaIndex  ]> 3* spatiere)                  { if there is a big gap between pointers then something must be wrong (one of the bases points to a wrong position way back into the chromatogram) }
 OR (Side.Ptr2Smpl[ctCromaIndex  ]> Side.Ptr2Smpl[ctCromaIndex+1]) then
  begin
   i:= Side.Ptr2Smpl[ctCromaIndex+1]- spatiere;
   Side.Ptr2Smpl[ctCromaIndex]:= i;
   Log.AddVerb('Invalid sample spaceation. Fixed._');
  end;
end;




function TAbiObj.CorrectBaseOrder(VAR Side: TSide): Boolean;                                     { Returns TRUE if invalid base pos found }
VAR NewPos, CurBase: Cardinal;
begin
 Result:= FALSE;

 { Cazul VIENA LAB }
 { Cazul in care primii doi pointeri sunt mai mari decat restul: 491 503, 257 267 281 298.       viennalab case/ bases shifted with 14\893_UTR-F1___.ab1  }
 { This must be FIRST! }
 for CurBase:= 4 downto 1 DO
  if Side.Ptr2Smpl[CurBase] > Side.Ptr2Smpl[CurBase+ 1]  then                                    { Pointer is shifted to right }
     begin
      Result:= TRUE;
      Side.Ptr2Smpl[CurBase]:= Side.Ptr2Smpl[CurBase+ 1] -1;                                     { Brute force solution }
      Log.AddVerb('Invalid pointer position at the beginning of the chromatogram (base no:' + IntToStr(CurBase)+ '). Fixed.');
     end;


 { Check current base agains the next base }
 for CurBase:= 1 to Side.NrOfBases-1 DO
  if Side.Ptr2Smpl[CurBase] > Side.Ptr2Smpl[CurBase+ 1] then                                     { Pointer is shifted to right }
   begin
    Result:= TRUE;
    Side.Ptr2Smpl[CurBase]:= Side.Ptr2Smpl[CurBase-1]+ spatiere;                                 { Brute force solution }
   end;

 { Also check the last base }
 if Side.Ptr2Smpl[Side.NrOfBases] < Side.Ptr2Smpl[Side.NrOfBases- 1] then
  begin
    NewPos:= Side.Ptr2Smpl[Side.NrOfBases- 1]+ spatiere;
    if NewPos > NoOfSamples
    then NewPos:= NoOfSamples;
    Side.Ptr2Smpl[Side.NrOfBases]:= NewPos;
    Result:= TRUE;
  end;

 if Result
 then Log.AddVerb('Invalid pointers position. Fixed.');     { This message appears to often so don't show it }
end;


procedure TAbiObj.CheckLastElem(VAR Side: TSide);
VAR correction: Cardinal;
begin
 if Side.Ptr2Smpl[Side.NrOfPtrs] <= Side.Ptr2Smpl[Side.NrOfPtrs- 1] then    { daca ultimul pointer e mai mic decat penultimul pointer }
  begin
    { Set the wrong pointer to the end of chroma }
    correction:= Side.Ptr2Smpl[Side.NrOfPtrs-1]+ spatiere;

    { Make sure I don't go out of chroma }
    if correction> NoOfSamples
    then correction:= NoOfSamples;

    Side.Ptr2Smpl[Side.NrOfPtrs]:= correction;
    Log.AddVerb('Last base (' + IntToStr(Side.NrOfPtrs) +'/'+ string(Side.BaseArray[Side.NrOfPtrs])+') was pointing to sample '+ IntToStr(Side.Ptr2Smpl[Side.NrOfPtrs]) +' instead of '+ IntToStr(correction) +'. Fixed.');
  end;
end;








{===============================================================================
                                LOG / DEBUGGING
===============================================================================}
function TAbiObj.ReturnAddresses: string;
begin
  {HEADER}
  Result:= '';
  Result:= Result
   + 'Header size: 28x '+ IntToStr(H.MainDirectory.NrOfElems)+ ' = '+ IntToStr(H.MainDirectory.NrOfElems* 28)+ 'Bytes'+ CRLF
   + 'Header count: '          + IntToStr(H.MainDirectory.NrOfElems)+ ' directories'+ CRLF
   + 'Header address: '        + IntToStr(H.MainDirectory.DataOffset)+ CRLF
   + 'Bases array ofset1: '    + IntToStr(Side1.BaseArrayOfset)+ CRLF
   + 'Bases array ofset2: '    + IntToStr(Side2.BaseArrayOfset)+ CRLF
   + 'Ptrs array offset1: '    + IntToStr(Side1.PtrsMXOffset)  + CRLF
   + 'Ptrs array offset2: '    + IntToStr(Side2.PtrsMXOffset)  + CRLF;

  {TACES}
  Result:= Result+ 'Trace order: '+ Trace1.Name+' '+ Trace2.Name+' '+ Trace3.Name+' '+ Trace4.Name+ LBRK;
  with Trace1 DO
   Result:= Result+ 'Trace '+ Name+ ' = '+IntToStr(NumberOf)+'samples. Address: '+ IntToStr(Offset)+ CRLF;
  with Trace2 DO
   Result:= Result+ 'Trace '+ Name+ ' = '+IntToStr(NumberOf)+'samples. Address: '+ IntToStr(Offset)+ CRLF;
  with Trace3 DO
   Result:= Result+ 'Trace '+ Name+ ' = '+IntToStr(NumberOf)+'samples. Address: '+ IntToStr(Offset)+ CRLF;
  with Trace4 DO
   Result:= Result+ 'Trace '+ Name+ ' = '+IntToStr(NumberOf)+'samples. Address: '+ IntToStr(Offset)+ CRLF;
end;



function TAbiObj.AbiProperties: string;
begin
 Result:=
     '  File version: '   + Tab + IntToStr(H.Version)       + CRLF
   + '  No of bases1: '   + Tab + IntToStr(Side1.NrOfBases) + CRLF
   + '  No of bases2: '   + Tab + IntToStr(Side2.NrOfBases) + CRLF
   + '  No of QVs 1: '    + Tab + IntToStr(Side1.QvCount)   + CRLF
   + '  No of QVs 2: '    + Tab + IntToStr(Side2.QvCount)   + CRLF
   + '  No of pointers1: '+ Tab + IntToStr(Side1.NrOfPtrs)  + CRLF
   + '  No of pointers2: '+ Tab + IntToStr(Side2.NrOfPtrs)  + CRLF
   + '  Tags: '           + Tab + string(H.AbiTags)+ CRLF;
end;



function TAbiObj.ShowDirectoriesRawData: string;
var i: Integer;
begin
 Result:= '';
 for i:= 0 TO H.MainDirectory.NrOfElems-1 DO
   Result:= Result+
     (string(Directories[i].Name)+ ':  '+
    IntToStr(Directories[i].ElemType) + '  '+
    IntToStr(Directories[i].ElemSize) + '  '+
    IntToStr(Directories[i].NrOfElems)+ '  '+
    IntToStr(Directories[i].DataSize) + '  '+
    IntToStr(Directories[i].DataOffset))+ CRLF;

 Result:= StringReplace(Result, #0, '?', [rfReplaceAll]);                                          { unele fisiere ABI corupte, baga caractere NULL ipe care trebuie sa le scot ca altfel nu pot sa le bag in log care e de tip TMemo si care nu accepta caract NULL }
end;



procedure TAbiObj.writeBases(Side: TSide; FileName: string);
VAR CurBase: Cardinal;
    s: string;
begin
 s:= '';
 for CurBase:= 1 to Side.NrOfBases DO
  s:= s+ string(Side.BaseArray[CurBase])+ ' '+ i2s(CurBase) + ' -> ' + i2s(Side.Ptr2Smpl[CurBase])+ CRLF;
 stringtofile(FileName, s, woOverwrite);
end;




end.
