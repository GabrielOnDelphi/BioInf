UNIT ScfBase;

{=============================================================================================================
 Gabriel Moraru
 2016.07
==============================================================================================================


 Base class for SCFRead/SCFWrite
 ===============================================================================

 Verified ok (except the SCF Writer)

 Documentation:
           Online documentatie for SCF version 3.0: http://staden.sourceforge.net/manual/formats_unix_8.html
           Current Code Set: IUPAC (NC-IUB). Code Set: (A,C,G,T,-)   (default)
           The 4 traces form together the Sample Matrix.
           uInt4 = Cardinal

 The content of the SCF MAY look like this:

           H.VERSION         3.00
           H.SampleSize      2
           H.spare[0 to 17]  RReadUInt4(FWriteOld);
           SamplesOffset     128
           NrOfSamples       8156
           BasesOffset       65376
           NrOfBases         632
           comments_offset   72960
           comments_size     326
           PrivateSize       0
           PrivateOffset     73286

 Swaping cardinals:
           Let's say we have 616 bases. The hex for it is $268.
           Unswapped it will be saved to disk as 68 02 00 00.
           Swapped it will be saved as 00 00 02 68.

===============================================================================}

INTERFACE

USES
   System.SysUtils,
   CubicDNA, ccStreamMem, ccCore, ccINIFile, clRamLog, ccRichLog;

TYPE
 { Packed record:  DON'T CHNAGE THE ORDER!
   The order of these fileds is critical when saving to disk. }
  RHeaderSCF= packed Record
    MagicNumber     : array[1..4] of AnsiChar;                                                     { should equal 779314022                    }{.scf}
    NoOfSamples     : Cardinal;                                                                    { Number of elements in Samples matrix      }{number of samples }
    SamplesOffset   : Cardinal;                                                                    { Byte offset from start of file            }{offset from the start of the file to the samples }
    NoOfBases       : Cardinal;                                                                    { Number of bases in Bases matrix           }
    ClipLeft        : Cardinal;                                                                    { OBSOLETE: No. bases in left clip (vector) }{bases left clip - we don't care }
    ClipRight       : Cardinal;                                                                    { OBSOLETE: No. bases in right clip (qual)  }{bases right clip - we don't care }
    BasesOffset     : Cardinal;                                                                    { Byte offset from start of file            }{bases offset - we don't care }
    CommentsSize    : Cardinal;                                                                    { Number of bytes in Comment section        }{comments size in bytes }
    CommentsOffset  : Cardinal;                                                                    { Byte offset from start of file            }{offset from the start of the file to the comments }
    Version         : array[1..4] of AnsiChar;                                                     { "version.revision", eg '3' '.' '0' '0'    }{3.00}
    SampleSize      : Cardinal;                                                                    { Size of samples in bytes 1 OR 2           }{8bit samples OR 16bit samples }
    Code_set        : Cardinal;                                                                    { code set used (but ignored!)              }{code_set - we don't care }
    PrivateSize     : Cardinal;                                                                    { No. of bytes of Private data, 0 if none   }{we don't care}
    PrivateOffset   : Cardinal;                                                                    { Byte offset from start of file            }{we don't care}
    Spare           : array[0..17] of Cardinal;                                                    { string[18]  Unused                        }{absolutely nothing}  {uInt4}
  end;

  
 TBase_v30= packed record                                                                          { Type definition for the sequence dat    }
    Ptr2Smpl : Cardinal;                                                                           { Pointer into Samples matrix for this base position}
    prob_A   : BYTE;                                                                               { Probability of it being an A            } {uInt1}
    prob_C   : BYTE;                                                                               { Probability of it being an C            }
    prob_G   : BYTE;                                                                               { Probability of it being an G            }
    prob_T   : BYTE;                                                                               { Probability of it being an T            }
    Base     : AnsiChar;                                                                           { Called base character                   }
    Spare    : array[0..2] of Byte;                                                                { Spare                                   }
  end;


 TBaseArray  = array of TBase_v30;
 TWordTrace  = array of Word;                                                                      { In SCF file, the samples could have negative value because there are stored as Delta values. Once the deltas are unpacked all values should be positive. }
 TDiskTrace  = array of Smallint;


 TScfObj = class(TObject)
  private
    procedure setNrOfBases  (Value: Integer);
    function  getNrOfBases: Integer;
    procedure setNrOfSamples(Value: Cardinal);
    function  getNrOfSamples: Cardinal;
  protected
    Log: TRamLog;
    FFileName: string;
    FStream: TCubicMemStream;
    procedure CheckQVs;
  public                                                                                           { SCF Structure: }
    H          : RHeaderSCF;
    TraceA     : TWordTrace;                                                                       { Indexed in 0 }{ 4 trace-rui formeaza un SAMPLE MATRIX. Pe 2 octeti desi in unele cazuri am nevoie doar de 1 byte. }
    TraceC     : TWordTrace;                                                                       { Indexed in 0 }{ 4 trace-rui formeaza un SAMPLE MATRIX }
    TraceG     : TWordTrace;                                                                       { Indexed in 0 }{ 4 trace-rui formeaza un SAMPLE MATRIX }
    TraceT     : TWordTrace;                                                                       { Indexed in 0 }{ 4 trace-rui formeaza un SAMPLE MATRIX }
    BaseArray  : TBaseArray;                                                                       { Indexed in 1 }{ MATRICEA BAZELOR, contine si QV MATRIX. Indexata in 1 }
    Comments   : string;
    PrivateData: string;

    constructor Create(aLog: TRamLog);                                                             { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
    destructor  Destroy; override;
    procedure Clear;

    function  Bases: BaseString;
    property  FileName   : string   read FFileName;
    property  NoOfBases  : Integer  read getNrOfBases   write setNrOfBases;
    property  NrOfSamples: Cardinal read getNrOfSamples write setNrOfSamples;
 end;




IMPLEMENTATION





{--------------------------------------------------------------------------------------------------
   CONSTRUCTOR
--------------------------------------------------------------------------------------------------}
constructor TScfObj.Create(aLog: TRamLog);                                                         { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
begin
 inherited Create;                                                                                 { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }
 Log:= aLog;
 FStream:= TCubicMemStream.Create;
end;

destructor TScfObj.Destroy;
begin
 FreeAndNil(FStream);
 inherited;
end;



procedure TScfObj.Clear;
begin
 FillChar(H, SizeOf(H), 0);                                                                        { In Delphi, FillChar fills Count contiguous bytes (referenced by X) with the value specified by Value (Value can be type Byte or Char). Warning:    This function does not perform any range checking.}
 FFileName   := '';
 Comments    := '';
 PrivateData := '';
 NoOfBases   := 0;
 NrOfSamples := 0;
 SetLength(BaseArray, 0);
 SetLength(TraceA, 0);
 SetLength(TraceC, 0);
 SetLength(TraceG, 0);
 SetLength(TraceT, 0);
end;





{--------------------------------------------------------------------------------------------------
   NrOfBases
--------------------------------------------------------------------------------------------------}
procedure TScfObj.setNrOfBases(Value: Integer);
VAR cNrOfBases: Cardinal;
begin
 cNrOfBases:= Value;
 if H.NoOfBases<> cNrOfBases then
  begin
   H.NoOfBases:= cNrOfBases;
   SetLength(BaseArray, H.NoOfBases+ IndexedIn1);                                                  { Indexed in 1 }
  end;
end;


function TScfObj.getNrOfBases: Integer;
begin
 Result:= H.NoOfBases;
end;




{--------------------------------------------------------------------------------------------------
   NrOfSamples
--------------------------------------------------------------------------------------------------}
procedure TScfObj.setNrOfSamples(Value: Cardinal);
VAR cardinalNrOfSamples: Cardinal;
begin
 cardinalNrOfSamples:= Value;                                                                      { Convert integer to Cardinal }
 if H.NoOfSamples<> cardinalNrOfSamples then
  begin
   H.NoOfSamples:= cardinalNrOfSamples;

   { Reserve memory for samples }
   SetLength(TraceA, H.NoOfSamples);                                                               { The matrix is indexed in zero }
   SetLength(TraceC, H.NoOfSamples);
   SetLength(TraceG, H.NoOfSamples);
   SetLength(TraceT, H.NoOfSamples);
  end;
end;


function TScfObj.getNrOfSamples: Cardinal;
begin
 Result:= H.NoOfSamples;
end;






{ Check against wrong input (Treat case Bradley 'PXY21-1_B_92294.scf' where QV values could be up to 255. }
procedure TScfObj.CheckQVs;
VAR
   BadQVs: Boolean;
   i: Integer;
begin
 BadQVs:= FALSE;
 for i:= IndexedIn1 to NoOfBases DO
  if (BaseArray[i].prob_A> 100) OR (BaseArray[i].prob_c> 100) OR (BaseArray[i].prob_g> 100) OR (BaseArray[i].prob_t> 100) then
   begin
    BadQVs:= TRUE;
    Log.AddWarn('Invalid QV info found (QV probability over 100%)! Attempting to recover data.');
    Break;
   end;

 if BadQVs then
  for i:= IndexedIn1 to NoOfBases DO
   begin
    BaseArray[i].prob_A:= BaseArray[i].prob_A DIV 5;
    BaseArray[i].prob_C:= BaseArray[i].prob_c DIV 5;
    BaseArray[i].prob_G:= BaseArray[i].prob_g DIV 5;
    BaseArray[i].prob_T:= BaseArray[i].prob_t DIV 5;
   end;
end;




function TScfObj.Bases: BaseString;
VAR i: Integer;
begin
 Result:= '';
 for i:= IndexedIn1 TO NoOfBases
  DO Result:= Result+ TBase(BaseArray[i].Base);
end;


end.
