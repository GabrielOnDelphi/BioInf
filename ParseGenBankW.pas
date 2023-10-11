UNIT ParseGenBankW;

{=============================================================================================================
 Gabriel Moraru
 2016.07
==============================================================================================================

   GBK to Fasta exporter
===============================================================================}

INTERFACE
USES
   Winapi.Windows, System.SysUtils, System.Classes,
   ccCore, ccINIFile, clRamLog, ccRichLog, ReadBasicSample, ParseGenBankR;

TYPE
 TGbkObjExp = class(TGbkObj)
   private
    function wrapOrigin: string;
   protected
   public
     Wrap: Boolean;                                                                               { If true, the bases are wrapped to 80 chars when saving to disk }
     NameUse_Strain  : Boolean;                                                                   { which fields should be included in file's name }
     NameUse_Clone   : Boolean;
     NameUse_Isolate : Boolean;
     NameUse_Gene    : Boolean;
     NameUse_Product : Boolean;
     NameUse_Organism: Boolean;
     NameUse_Version : Boolean;
     NameUse_Definit : Boolean;
     IgnoreProduct   : Boolean;                                                                   { When generating filename, ignore the 'PRODUCT' field but only if 'GENE' field is not empty }
     constructor Create(aLog: TRamLog);
     destructor Destroy; override;
     procedure  Clear;   override;
     {Export}
     function  AsText: string;                                                                    { Get the content of this GBK object as text so I can save it to disk }
     procedure SaveAsFasta (CONST FullFileName: string);
     procedure SaveAsFasta_AutoName (CONST Folder: string);
     procedure GenerateComment;  override;
     function  GenerateNameFromField: string;
     function  GenerateNameFromFile (CONST Nume: string): string;
 end;





IMPLEMENTATION

USES
   ReadFasta, ccIO;







{--------------------------------------------------------------------------------------------------
   CREATE / DESTROY
--------------------------------------------------------------------------------------------------}
constructor TGbkObjExp.Create;
begin
 inherited Create(aLog);                                                                           { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }
 NameUse_Strain  := TRUE;                                                                          { which fields should be included in file's name }
 NameUse_Clone   := TRUE;
 NameUse_Isolate := TRUE;
 NameUse_Gene    := TRUE;
 NameUse_Product := TRUE;
 NameUse_Organism:= TRUE;
 NameUse_Version := TRUE;
 NameUse_Definit := TRUE;
 IgnoreProduct   := TRUE;                                                                          { When generating FileNameParent, ignore the 'PRODUCT' field but only if 'GENE' field is not empty }
 Wrap:= TRUE;
end;



destructor TGbkObjExp.Destroy;
begin
 inherited Destroy;
end;


procedure TGbkObjExp.Clear;
begin
 inherited;
 NameUse_Strain  := TRUE;                                                                          { which fields should be included in file's name }
 NameUse_Clone   := TRUE;
 NameUse_Isolate := TRUE;
 NameUse_Gene    := TRUE;
 NameUse_Product := TRUE;
 NameUse_Organism:= TRUE;
 NameUse_Version := TRUE;
 NameUse_Definit := TRUE;
 IgnoreProduct   := TRUE;                                                                          { When generating FileNameParent, ignore the 'PRODUCT' field but only if 'GENE' field is not empty }
end;







{--------------------------------------------------------------------------------------------------
   EXPORT AS GBK
--------------------------------------------------------------------------------------------------}

function TGbkObjExp.AsText: string;                                                                 { Get the content of this GBK object as text so I can save it to disk }
VAR
   i: Integer;

 procedure AddRecord(RecordName, Data: string; Force: Boolean = FALSE);                             { Force = force title to be written to disk even if Line is empty }
 VAR
    i: Integer;
    TSL: TStringList;
 begin
  RecordName:= MakeStringLongRight(RecordName, ' ', HeaderLength);

  if Length(Data) > 68 then
   begin
    Data:= WrapText(Data, 68);
    TSL:= SplitText(Data, CRLF);
    for i:= 1 to tsl.Count-1 DO
     tsl[i]:= '            '+ tsl[i];      // 12 spaces
    Data:= TSL.Text;
    FreeAndNil(TSL);
   end;
  Data:= RemoveLastEnter(Data);

  if Force
  OR (Data > '')
  then Result:= Result+ RecordName+ Data+ CRLF;
 end;

 procedure AddFeatField(FeatName, Data: string);                                 { Force = force title to be written to disk even if Line is empty }
 begin
  if Data > '' then
   begin
    //FeatName:= '/'+ FeatName;
    FeatName:= StringOfChar(' ', DataStartFeatures) + '/'+ FeatName;
    Result:= Result+ FeatName+ '="'+ Data+ '"' + CRLF;
   end;
 end;

begin
 Result:= '';

 AddRecord('LOCUS'          ,  Records.LOCUS);
 AddRecord('DEFINITION'     ,  Records.DEFINITION);
 AddRecord('ACCESSION'      ,  Records.ACCESSION);
 AddRecord('VERSION'        ,  Records.VERSION);
 AddRecord('KEYWORDS'       ,  Records.KEYWORDS+ '.');
 AddRecord('SOURCE'         ,  Records.Source.default, Records.Source.organism.default > '');
 AddRecord('  ORGANISM'     ,  Records.Source.organism.default);
 AddRecord('            '   ,  Records.Source.organism.Data+ '.');

 AddRecord('REFERENCE',  Records.REFERENCE.Default);
  for i:= 0 to High(Records.REFERENCE.RefList) DO
   begin
    AddRecord('  AUTHORS'    ,  Records.REFERENCE.RefList[i].authors);
    AddRecord('  TITLE'      ,  Records.REFERENCE.RefList[i].TITLE);
    AddRecord('  JOURNAL'    ,  Records.REFERENCE.RefList[i].JOURNAL);
    AddRecord('  PUBMED'     ,  Records.REFERENCE.RefList[i].PUBMED);
   end;

 AddRecord('COMMENT'        ,  Records.COMMENT);

 AddRecord('FEATURES'       ,  Records.FEATURES.Default);
 AddFeatField('organism'    ,  Records.FEATURES.organism);
 AddFeatField('isolate'     ,  Records.FEATURES.isolate);
 AddFeatField('clone'       ,  Records.FEATURES.clone);
 AddFeatField('mol_type'    ,  Records.FEATURES.mol_type);
 AddFeatField('strain'      ,  Records.FEATURES.strain);
 AddFeatField('note'        ,  Records.FEATURES.note);
 AddFeatField('gene'        ,  Records.FEATURES.gene);
 AddFeatField('codon_start' ,  Records.FEATURES.codon_start);
 AddFeatField('transl_table',  Records.FEATURES.transl_table);
 AddFeatField('product'     ,  Records.FEATURES.product);
 AddFeatField('protein_id'  ,  Records.FEATURES.protein_id);
 AddFeatField('db_xref'     ,  Records.FEATURES.db_xref);
 AddFeatField('translation' ,  Records.FEATURES.translation);

 AddRecord('ORIGIN'         ,  '', TRUE);                     { This is the sequence }
 Result:= Result+ WrapOrigin;  //string(Records.ORIGIN);
 Result:= Result+ ctGbkEndOfSample;
end;




function TGbkObjExp.WrapOrigin: string;
VAR
   Line, s: string;
   i: Integer;

   procedure AddLine;
   VAR LineHdr: string;
   begin
    if Length(Line) = 60
    then LineHdr:= IntToStr(i- 60 +1)+ ' '
    else LineHdr:= IntToStr(i- Length(Line))+ ' ';

    LineHdr:= MakeStringLongLeft(LineHdr, ' ', 10);

    Line:= InsertCharEvery(' ', Line, 10);
    Result:= Result+ LineHdr+ Line+ CRLF;                                     { cand textul are fix 'spnSplit.Value' caractere pe rand, imi adauga un rand gol. }
    Line:= '';
   end;

Begin
 Result:= '';
 Line:= '';
 s:= String(Records.ORIGIN);

 for i:= 1 TO Length(s) DO
  begin
   Line:= Line+ s[i];                                           {TODO 5: THIS IS SLOW. I need some kind of string builder here }
   {
   if i mod 10 = 0
   then Line:= Line+ ' '; }

   if i mod 60 = 0       { trunchez liniile la 65 caractere }
   then AddLine;
  end;

  if FirstChar(Line) <> ' '
  then AddLine;
End;








{--------------------------------------------------------------------------------------------------
   EXPORT AS FASTA
--------------------------------------------------------------------------------------------------}
procedure TGbkObjExp.SaveAsFasta (CONST FullFileName: string);
VAR ObjFasta: TFastaObj;
begin
  ObjFasta:= TFastaObj.Create(RamLog);
  TRY
   ObjFasta.Bases:= Bases;
   ObjFasta.Comment:= Comment;

   ObjFasta.SaveAsFasta(FullFileName, Wrap);
  FINALLY
   FreeAndNil(ObjFasta);
  END; 
end;


procedure TGbkObjExp.SaveAsFasta_AutoName(CONST Folder: string);
VAR ObjFasta: TFastaObj;
    s: String;
begin
 ObjFasta:= TFastaObj.Create(RamLog);
 TRY
  s:= GenerateNameFromField;
  FileName:= Trail(Folder)+ s+ '.FASTA';

  ObjFasta.Bases:= Bases;
  ObjFasta.Comment:= s;

  { Check name lenght }
  if Length(FileName) < Max_Path
  then ObjFasta.SaveAsFasta(FileName, Wrap)                                          { TRY TO SAVE AS FASTA }
  else RamLog.AddWarn(' File name is too long: '+ FileName)
 FINALLY
  FreeAndNil(ObjFasta);
 END;
end;













{--------------------------------------------------------------------------------------------------
   GENEREAZA NUMELE
--------------------------------------------------------------------------------------------------}

function TGbkObjExp.GenerateNameFromField: string;                                   { cand am multi GenBank }
VAR
   strain, definit, isolate,
   clone, gene, product,
   organism, version: string;
begin
 if (Records.FEATURES.organism <> '') AND NameUse_Organism                           { which Records should be included in file's name }
 then organism:= Records.FEATURES.organism;

 if (Records.DEFINITION <> '')
 AND NameUse_Definit
 then definit := '  ' + Records.DEFINITION;

 if (Records.FEATURES.strain <> '')
 AND NameUse_Strain
 then strain := '  strain ' + Records.FEATURES.strain;

 if (Records.FEATURES.clone <> '')
 AND NameUse_Clone
 then clone := '  clone ' + Records.FEATURES.clone;

 if (Records.FEATURES.isolate <> '')
 AND NameUse_Isolate
 then isolate := '  isolate '+ Records.FEATURES.isolate;

 if (Records.FEATURES.gene    <> '')
 AND NameUse_Gene
 then gene:= '  '+ Records.FEATURES.gene;

 if (Records.FEATURES.product <> '')
 AND NameUse_Product
 then
   if (Records.FEATURES.gene <> '')
   AND IgnoreProduct AND NameUse_Gene                              { When generating FileName, ignore the 'PRODUCT' field but only if 'GENE' field is not empty }
   then product := ''
   else product := '  ' + Records.FEATURES.product;

 if (Records.VERSION <> '')
 AND NameUse_Version
 then version := '  ' + Records.VERSION;

 Result:= organism+ definit+ strain+ clone+ isolate+ gene+ product+ version;
 Result:= CorrectFilename(Result, ' ');
end;



function TGbkObjExp.GenerateNameFromFile(CONST Nume: string): string;                                 { mai imi trebuie asta? poate pt Cris }
VAR
  strain, isolate, clone: string;
begin
 if Records.FEATURES.strain <> '' then strain := ' strain ' + Records.FEATURES.strain;
 if Records.FEATURES.clone  <> '' then clone  := ' clone '  + Records.FEATURES.clone;
 if Records.FEATURES.isolate<> '' then isolate:= ' isolate '+ Records.FEATURES.isolate;

 Result:= ExtractOnlyName(Nume);
 Result:= CopyTo(Result,1,Pos(' ',Result));
 Result:= Records.FEATURES.organism+ strain+ clone+ isolate+ ' [] '+ Result + Records.VERSION;
 Result:= CorrectFilename(Result, ' ');
end;



procedure TGbkObjExp.GenerateComment;           { GENERATE COMMENTS }
begin
 Comment:= '> '+ GenerateNameFromField;
end;









end.
