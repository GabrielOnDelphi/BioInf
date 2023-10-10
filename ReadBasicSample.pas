
UNIT ReadBasicSample;

{===========================================================================================
 Heracle BioSoft SRL
 2016.04.28
 Common ancestor for: TCubeAbstract, TFastaObj, TGbkObj
============================================================================================}

INTERFACE

USES System.SysUtils, System.StrUtils, ccCore, ccINIFile, CubicDNA, clRamLog, ccRichLog;

CONST
   ctDefaultRowLengt= 70;                                                                  { Default row length for FASTA file. After this length the BASES string will be broken to a new row }

TYPE
   RDetectedVector= record
      Starts: Integer;                                                                     { The right vector will start at position 10  (for example) and end at 50 }
      Ends  : Integer;                                                                     { The right vector will start at position 100 (for example) and end at 90 (backwords) }
      Name: string;
      Value: BaseString;
   end;

   RVectorData = record                                                                    { TCubeOnj and TSimpleSeq keeps the info about the vectors in this record }
     Detector : TObject;
     Left     : RDetectedVector;
     Right    : RDetectedVector;
    end;

TYPE
  TBasicSample = class(TObject)
   private
     FParentType: TBioFileType;
     FComment: String;                                                                     { Pe HDD fisierele SEQ sunt stocate ca o cloectie de siruri separate de CRLF }
     procedure setComment (CONST Value: string);
     procedure setFileName(const Value: string);
   protected
     RamLog: TRamLog;
     FBases: BaseString;
     function  getFileName: string;               virtual;                                 { This will be overriden in TCubeObj }
     function  getNoBases: Integer;               virtual;
     procedure setNoBases (CONST Value: Integer); virtual;
     procedure GenerateComment;                   virtual;
   public
     Vectors     : RVectorData;
     IsPart      : Boolean;                                                                { True if this cube is a FASTA which is part of a multi-FASTA object }
     CleanedVectors: string;                                                               { It will not be created implicitelly by 'Create' but it will be freed by 'Destroy' }
     FFileName   : string;
     constructor Create(aLog: TRamLog);
     destructor Destroy; override;
     procedure  Clear; virtual;
     function   ScreenName: string;
     function   ShortName: string;
     function   CommentIsEmpty: Boolean;

     procedure  RemoveAllGaps;
     function   DetectedVectors: string;                                                   { Return a string containing the name of the Left and Right vectors }
     procedure  VectorsClearColor; virtual;                                                { All existent cells already marked as Vectors, are cleared. List of detected vectors is also cleared }
     procedure  DetectVectors;  virtual;                                                   { Detect vectors. Set cell colors (blue) for the detected vectors. }
     function   Search (CONST SearchStr: BaseString; Offset: Integer): Integer; virtual;

     property   Comment        : String       Read FComment      Write setComment;
     property   ParentType     : TBioFileType Read FParentType   Write FParentType;        { What is the original provenience (input) of this Cube. Example: SCF, ABI, SEQ, FASTA, etc. Needed if the user wants to save this seq back to file (after he edited the seq).    The FParentType is the name of the file from which this object was loaded. Make sense only if the parent is a multi file (multiFasta, multiGBK) }
     property   FileName       : string       Read getFileName     Write setFileName;      { The actual name. The user may change it. For contigs, it only contains the filename, without path. Path is extracted from ProjectName }
     property   NoOfBases      : Integer      Read getNoBases    Write setNoBases;         { are grija sa adune un 1 pt ca CellMX e indexat in 1 nu in 0 }  { Number of bases in Bases matrix.  IMPORTANT. BasesNr e diferit de BasesNrOrig pentru ca userul/asamblarea poate adauga baze noi in plus. Cand import un fisier o sa folosesc BasesNrOrig, insa dupa aceea o sa fol. BasesNr }
  end;


IMPLEMENTATION

USES
   VectorProcessing, ccIO;








{--------------------------------------------------------------------------------------------------
   CONSTRUCTOR
--------------------------------------------------------------------------------------------------}
constructor TBasicSample.Create(aLog: TRamLog);                                                           { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
begin
 inherited Create;                                                                                        { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }
 Vectors.Detector:= NIL;
 RamLog:= aLog;
 Assert(aLog<> NIL, 'TBasicSample.Log is NIL!');

 FParentType:= bfNone; { This must stay here }

 { VECTORS }
 VectorsClearColor;
end;


destructor TBasicSample.Destroy;
begin
 inherited;
end;


procedure TBasicSample.Clear;
begin
 CleanedVectors  := '';
 FBases          := '';
 FComment        := '';
 FParentType     := bfNone;
 IsPart          := FALSE;                                                                         { Shows if this GBK is part of a multi-GBK object }
//DELETE! DEL DEL DEL RamLog.Clear;  { This will not clear RichLog. The RichLog might have text from other RamLogs. The the GUI control the RichLog }

 { VECTORS }
 VectorsClearColor;
end;




procedure TBasicSample.setFileName(const Value: string);
begin
 FFileName := Value;
 ParentType:= GetParentObj(Value);
end;


function TBasicSample.getFileName: string;
begin
 Result:= FFileName;
end;




procedure TBasicSample.setNoBases(CONST Value: Integer);
begin

 raise Exception.Create('TBasicSample.setNoBases is abstract!');
end;



function TBasicSample.getNoBases: Integer;
begin
 Result:= Length(FBases);
end;






{--------------------------------------------------------------------------------------------------
   COMMENTS
--------------------------------------------------------------------------------------------------}
procedure TBasicSample.GenerateComment;
begin
 FComment:= '> '+ ScreenName;
end;



procedure TBasicSample.setComment(CONST Value: String);
begin
 FComment:= ValidateComment(Value);                                                                { Replace enter with a space, Scoate caracterele ascii sub 32, Force '>' at the begining of the string }
end;



function TBasicSample.CommentIsEmpty: Boolean;
VAR s: string;
begin
 s:= FComment;
 s:= ReplaceString(s, ' ', '');
 s:= ReplaceString(s, '>', '');
 s:= TrimEnters(s);
 Result:= s= '';
end;





{--------------------------------------------------------------------------------------------------
                                     STUFF
--------------------------------------------------------------------------------------------------}
function TBasicSample.ScreenName: string;
begin
 Result:= ExtractOnlyName(FileName)
end;



function TBasicSample.ShortName: string;
begin
 Result:= ExtractFileName(FileName);
end;



function TBasicSample.Search(CONST SearchStr: BaseString; Offset: Integer): integer;
begin
 Result:= PosEx(String(SearchStr), String(FBases), Offset);
 {TODO: Search: IMPROVE SPEED HERE! }
end;



procedure TBasicSample.RemoveAllGaps;
begin
 FBases:= ReplaceString(FBases, GAP, '');
end;






{--------------------------------------------------------------------------------------------------
   VECTORS
--------------------------------------------------------------------------------------------------}

{ NOTE:
  Vectors are detected everytime the cub is edited/changed by calling 'DetectVectors' in 'getBases' function. }
procedure TBasicSample.VectorsClearColor;
begin
 Vectors.Right.Starts:= -1;
 Vectors.Right.Ends  := -1;
 Vectors.Right.Name  := '';
 Vectors.Right.Value := '';
 Vectors.Left.Starts := -1;
 Vectors.Left.Ends   := -1;
 Vectors.Left.Name   := '';
 Vectors.Left.Value  := '';
end;


procedure TBasicSample.DetectVectors;                                                              { WEIRD STUFF TO REMEMBER: If the vector contains an IUPAC base it will work (the vector will be detected and cut) because we extrapolate the vectors. However, if the CONTIG contains an IUPAC base in the vector region, the vector won't be detected because... well... there is no match between contig and the vector. That IUPAC base breaks the match. }
begin
 if Vectors.Detector<> NIL
 then (Vectors.Detector as TVectorDetector).DetectVectors(Self);
end;


function TBasicSample.DetectedVectors: string;                                                     { Return a string containing the name of the Left and Right vectors }
begin
 if Vectors.Detector= NIL then EXIT('');

 if Vectors.Left.Starts> 0
 then Result:= Vectors.Left.Name;

 if Vectors.Right.Starts> 0
 then
   if Result= ''
   then Result:= Vectors.Right.Name
   else Result:= Result+ ', '+ Vectors.Right.Name;
end;



END.
