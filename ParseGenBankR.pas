
UNIT ParseGenBankR;

{===============================================================================
 Heracle BioSoft
 2016.04.13
 
 Documentation GBK format:
    http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
    http://www.nmpdr.org/FIG/wiki/view.cgi/FIG/GBK
 ==============================================================================}

INTERFACE
USES
   System.StrUtils, System.SysUtils, System.Classes,
   ccCore, ccINIFile, CubicDNA, ReadBasicSample, clRamLog, ccRichLog;

CONST
   HeaderLength= 12;                                                                              { the content of every filed starts at column 12+1 }
   ctGbkEndOfSample= '//';
   DataStart= 13;
   DataStartFeatures= 22;
   DataStartBases= 11;

TYPE
  TFeatures= record
    Default     : string;
    organism    : string;
    isolate     : string;
    clone       : string;
    mol_type    : string;
    strain      : string;
    note        : string;
    gene        : string;                                                                         //  /gene="pmoA"
    codon_start : string;                                                                         //  /codon_start=1
    transl_table: string;                                                                         //  /transl_table=11
    product     : string;                                                                         //  /product="particulate methane monooxygenase protein A"
    protein_id  : string;                                                                         //  /protein_id="AAF08213.1"
    db_xref     : string;                                                                         //  /db_xref="GI:6424934"
    translation : string;                                                                         //  Poate fi pe mai multe randuri:   /translation="GDWDFWVDWKDRRMWPTVVPILGVTFAAASQAFWWVNFRLPFGAVFAALGLLIGEWINRYVNFWGWTYFPISLVFPSALIVPAIWLDVILLLSGSYVITAVVGSLGWGLLFYPNNWPAIAAFHQATEQHGQLMTLADLIGFHFVRTSMPEYIRMVERGTL
  end;

  ROrganism= record
   default: string;
   Data: string;
  end;

  TSource= record
    Default     : string;
    Organism    : ROrganism;
  end;

  TReference= record                                                                              { GBK can have more than one reference }
    authors      : string;
    title        : string;
    journal      : string;
    pubmed       : string;
  end;

  RReferences= record
    Default  : string;
    RefList  : array of TReference;                                                               { GBK can have more than one reference }
  end;

  TBank= record
    LOCUS        : string;
    DEFINITION   : string;
    ACCESSION    : string;
    VERSION      : string;
    KEYWORDS     : string;
    ORIGIN       : BaseString;                                                                    { Aici e secventa propriu zisa }
    COMMENT      : string;
    SOURCE       : TSource;
    FEATURES     : TFeatures;
    REFERENCE    : RReferences;
  end;

 TGbkObj = class(TBasicSample)
   private
     FBody: TStringList;                                                                          { This text is temporary. Parse will irreversibly clear the text. }
     function  ExtractRecord   (RecordName: string): TStringArray;
     function  ParseSimpleRec  (RecordName: string): string;
     function  ExtractDataFrom (CONST Line: string; StartPos: Integer= DataStart): string;
     function  FindField       (StrArray: TStringArray; FieldName: string): Integer;
   protected
     function  getNoBases: Integer; override;
     function  ParseFeaturesSub (StrArray: TStringArray; CONST s: string): string;
     procedure ParseFeatures;
     procedure ParseSource;                                                                       { Parse SOURCE }
     procedure ParseOrigin;
   public
     Records: TBank;
     constructor Create(aLog: TRamLog);
     destructor Destroy; override;
     procedure  Clear;   override;
     procedure AddRawData(sLine: string);
     procedure Parse;

     function  Bases: BaseString;
     function  LoadFromFile (FullFileName: string): string;
 end;


function BreakSourceToSubfields(SourceField: string): TStringList;                                { Takes the "Source.Organism" field and breaks it in subfields. Places each subfield on a new line in a TSL. DON'T FORGET TO FREE THE OBJECT RETURNED AFTER YOU FINISH WITH IT. Cine il foloseste? PolyPro? }

IMPLEMENTATION

USES
   ccIO;

CONST
   RecordStart= 1;
   FieldStart = 3;                                                                                { This is the start of sub-items in a 'normal' record }




{--------------------------------------------------------------------------------------------------
   CREATE / DESTROY
--------------------------------------------------------------------------------------------------}
constructor TGbkObj.Create;
begin
 inherited Create(aLog);                                                                                      { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }
 FBody:= TStringList.Create;
 Clear;                                                                                          { Asta e intotdeauna ultimul }
end;


destructor TGbkObj.Destroy;
begin
 FreeAndNil(FBody);
 inherited Destroy;
end;


procedure TGbkObj.Clear;
begin
 inherited;
 FBody.Clear;

 WITH Records DO
  begin
   LOCUS           := '';
   DEFINITION      := '';
   ACCESSION       := '';
   VERSION         := '';
   KEYWORDS        := '';
   ORIGIN          := '';                                                                        { Aici e secventa propriu zisa:  1 ggggactggg acttctgggt tgactggaag gatcgtcgta tgtggccgac ggtcgtgccg }

   with SOURCE DO
    begin
     default         := '';
     organism.default:= '';
     organism.Data   := '';
    end;

   SetLength(REFERENCE.RefList, 0);
   REFERENCE.default:= '';                                                                       { REFERENCE }

   with FEATURES DO
    begin
     organism       := '';
     isolate        := '';
     clone          := '';
     mol_type       := '';
     strain         := '';
     note           := '';
     gene           := '';                                                                       //  /gene="pmoA"
     codon_start    := '';                                                                       //  /codon_start=1
     transl_table   := '';                                                                       //  /transl_table=11
     product        := '';                                                                       //  /product="particulate methane monooxygenase protein A"
     protein_id     := '';                                                                       //  /protein_id="AAF08213.1"
     db_xref        := '';                                                                       //  /db_xref="GI:6424934"
     translation    := '';                                                                       //   Poate fi pe mai multe randuri:   /translation="GDWDFWVDWKDRRMWPTVVPILGVTFAAASQAFWWVNFRLPFGAVFAALGLLIGEWINRYVNFWGWTYFPISLVFPSALIVPAIWLDVILLLSGSYVITAVVGSLGWGLLFYPNNWPAIAAFHQATEQHGQLMTLADLIGFHFVRTSMPEYIRMVERGTL
    end;
  END;
end;







{--------------------------------------------------------------------------------------------------
   PARSE
--------------------------------------------------------------------------------------------------}

procedure TGbkObj.Parse;                                                                          { Process raw data and load info into the appropriate Records }
begin
 Assert(FBody.Text > '', 'FBody.Text is empty!');

 Records.LOCUS      := ParseSimpleRec('LOCUS');
 Records.DEFINITION := ParseSimpleRec('DEFINITION');
 Records.ACCESSION  := ParseSimpleRec('ACCESSION');
 Records.VERSION    := ParseSimpleRec('VERSION');
 Records.KEYWORDS   := ParseSimpleRec('KEYWORDS');
 Records.KEYWORDS   := RemoveLastChar(Records.KEYWORDS, '.');
 ParseSource;
 Records.REFERENCE.Default := ParseSimpleRec('REFERENCE');
 ParseFeatures;
 ParseOrigin;

 FBody.Clear;
end;






{--------------------------------------------------------------------------------------------------
   RECORD/FIELD PARSERS
--------------------------------------------------------------------------------------------------}

function TGbkObj.ExtractRecord(RecordName: string): TStringArray;                                     { Returns all lines belonging to this field }
VAR
   i, j: Integer;
   s: string;

  procedure AddElement(NewLine: string);
  begin
    SetLength(Result, Length(Result)+1);
    Result[High(Result)]:= NewLine;
  end;

begin
 SetLength(Result, 0);

 for i:= 0 to FBody.Count-1 DO
  begin
   s:= FBody[i];
   if Pos(RecordName, s)= 1 then                                                                      { Main recors always start at position 1. Subrecords always start at post 3 }
    begin
     AddElement(s);

     { Read next lines }
     if (i+1 = FBody.Count) then EXIT;                                                                { Don't go out of FBody }

     for j:= i+1 to FBody.Count-1 DO
       if FirstChar(FBody[j]) = ' '                                                                   { Subrecords always start at post 3 (start with 2 spaces) }
       then AddElement(FBody[j])
       else EXIT;
    end;
  end;
end;



function TGbkObj.FindField(StrArray: TStringArray; FieldName: string): Integer;                       { Returs the position (row) where the field was found }
VAR
   j: Integer;
   s: string;
begin
 Result:= -1;
 for j:= 0 to High(StrArray) DO
  begin
   s:= StrArray[j];
   s:= Copyto(s, FieldStart, HeaderLength);
   s:= RemoveSpaces(s);
   if s= FieldName
   then EXIT(j);
  end;
end;



function TGbkObj.ParseSimpleRec(RecordName: string): string;                                          { Returns the data contained in this record. This is for records that have no sub-fields! }
VAR
   i: Integer;
   StrArray: TStringArray;
begin
 Result:= '';
 StrArray:= ExtractRecord(RecordName);
 for i:= 0 to High(StrArray)
   DO Result:= Result+ extractDataFrom(StrArray[i])+ ' ';

 Result:= RemoveLastSpace(Result);
end;



function TGbkObj.extractDataFrom(CONST Line: string; StartPos: Integer= DataStart): string;           { Extract data from a row }
begin
 Result:= system.COPY(Line, StartPos, MAXINT);
end;











{-----------------
   Parse ORIGIN
------------------}
procedure TGbkObj.ParseOrigin;
VAR
   i: Integer;
   StrArray: TStringArray;
   Bases: string;
begin
 StrArray:= ExtractRecord('ORIGIN');
 for i:= 0 to High(StrArray)
   DO Bases:= Bases+ ExtractDataFrom(StrArray[i], DataStartBases);

 Bases:= RemoveLastSpace(Bases);
 Records.ORIGIN:= CubicDNA.CleanSequenceFormatings(BaseString(Bases));
end;




{-----------------
   Parse SOURCE
------------------}
procedure TGbkObj.ParseSource;
VAR
   Poz, i: Integer;
   StrArray: TStringArray;
   Line: string;
begin
 StrArray:= ExtractRecord('SOURCE');
 if Length(StrArray) = 0 then EXIT;

 Records.Source.default:= ExtractDataFrom(StrArray[0]);
 Poz:= FindField(StrArray, 'ORGANISM');

 if Poz> -1 then                                   { Read next row and see if it is an extension of the current row }
  begin
     Line:= ExtractDataFrom(StrArray[Poz]);
     Records.Source.organism.Default:= RemoveLastChar(Line, '.');

     Line:= '';
     if Poz+ 1 <= High(StrArray) then
      begin
       for i:= Poz+1 to High(StrArray)
        DO Line:= Line+ ExtractDataFrom(StrArray[i])+ ' ';

       Line:= RemoveLastSpace(Line);
       Line:= RemoveLastChar (Line, '.');
       Records.Source.organism.Data:= Line;                  { Scot punctul de la sfarsitul randului }   { Atentie, acest punct s-ar putea sa nu apara in fisierele generate de ARB }
      end;
 end;
end;


function BreakSourceToSubfields(SourceField: string): TStringList; { Takes the "Source.Organism" field which may look like this "Bacteria; Proteobacteria; environmental samples" and breaks it in subfields. Places each subfield on a new line in a TSL }
VAR s: string;
    CurPos: Integer;
begin
 Result:= TStringList.Create;
 for CurPos:= 1 to Length(SourceField)+1 DO                                                        { +1 ca sa ii dau sansa sa scrie si ultimul sub-field care NU se termina cu ';' }
  if  (CurPos<= Length(SourceField))
  AND (SourceField[CurPos]<> ';')
  then s:= s+ SourceField[CurPos]
  else
    begin
     Result.Add(Trim(s));
     s:= '';
    end;
end;





{---------------------
   Parse FEATURES
----------------------}
procedure TGbkObj.ParseFeatures;
VAR
   StrArray: TStringArray;
begin
 StrArray:= ExtractRecord('FEATURES');
 if Length(StrArray) > 0 then      { Get sub-FEATURES }
  begin
   Records.FEATURES.default := ExtractDataFrom(StrArray[0], DataStartFeatures);

   Records.FEATURES.organism:= ParseFeaturesSub (StrArray, '/organism=');
   Records.FEATURES.strain  := ParseFeaturesSub (StrArray, '/strain=');
   Records.FEATURES.isolate := ParseFeaturesSub (StrArray, '/isolate=');
   Records.FEATURES.clone   := ParseFeaturesSub (StrArray, '/clone=');
   Records.FEATURES.gene    := ParseFeaturesSub (StrArray, '/gene=');
   Records.FEATURES.product := ParseFeaturesSub (StrArray, '/product=');
  end;
end;


function TGbkObj.ParseFeaturesSub (StrArray: TStringArray; CONST s: string): string;
VAR i: Integer;
    MarkBegin, MarkEnd: Integer;
begin
 Result:= '';
 for i:= 0 TO High(StrArray) DO
  if Pos(s,StrArray[i]) >0 then
   begin
    MarkBegin:= Pos  ('"', StrArray[i]);
    MarkEnd  := PosEx('"', StrArray[i], MarkBegin+1);

    if MarkEnd= 0 then                                                                              { daca MarkEnd=0 inseamna ca... }
     begin
      MarkEnd:= PosEx('"', StrArray[i+1], 1);
      if MarkEnd=0
      then RamLog.AddError('[ERROR in ParseFeaturesSub]  Can''t find field delimiter for '+ s)
      else
       begin
        Result:= CopyTo(StrArray[i], MarkBegin+1, High(Integer))+ CopyTo(StrArray[i+1], 1, MarkEnd-1);        { ...acest camp e rupt pe doua sau mai multe randuri }
        Break;
       end;
     end
    else
     begin
      Result:= CopyTo(StrArray[i], MarkBegin+1, MarkEnd-1);
      break;
     end;                                                                                           { ...acest camp nu e rupt pe doua randuri }
   end;
end;



















{--------------------------------------------------------------------------------------------------
   I/O
--------------------------------------------------------------------------------------------------}
function TGbkObj.LoadFromFile (FullFileName: string): string;
begin
 Clear;
 FileName:= FullFileName;
 FBody.Text:= String(StringFromFileA(FullFileName));
 Parse;                                                                                            { Process raw data and load info into the appropriate Records }
end;


procedure TGbkObj.AddRawData(sLine: string);                                                       { Adds a line of data. The ReadGbkMulti will add multiple lines }
begin
 FBody.Add(sLine);
end;






{--------------------------------------------------------------------------------------------------
   OVERRIDES
--------------------------------------------------------------------------------------------------}
function TGbkObj.getNoBases: Integer;
begin
 Result:= Length(Records.Origin);
end;


function TGbkObj.Bases: BaseString;
begin
 Result:= Records.Origin;
end;








end.
