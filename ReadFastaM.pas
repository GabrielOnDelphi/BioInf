
UNIT ReadFastaM;

{=======================================================================================================================
 Heracle BioSoft SRL
 2016.07.13

=======================================================================================================================}

INTERFACE
USES
   Winapi.Windows, System.Classes, System.Types, System.SysUtils, System.Contnrs,
   ccCore, ccINIFile, CubicDNA, clRamLog, ccRichLog, ReadFasta, ReadBasicSample;

TYPE
 TMultiFasta = class(TObject)
  private
    function getFasta(Index: Integer): TFastaObj;
  protected
    List             : TObjectList;
    RamLog           : TRamLog;
    procedure generateVirtualNames;                                                                { When I load a MultiFasta file, all cubs will have the same name, so here I generate a unique name for each cub }
  public
    FileName         : string;
    CleanedVectors   : string;                                                                     { It will not be created implicitelly by 'Create' but it will be freed by 'Destroy' }
    RowLength        : Integer;
    Name_CustomName  : string;
    Name_UseIncrement: boolean;                                                                    { applies only if the input file contains more than one subsamples }
    Name_UseOrigName : boolean;
    Name_UseComments : boolean;
    OutputFolder     : string;                                                                     { The folder where the output files will be saved when the Split2Fasta function is used }
    AllowProtein     : Boolean;                                                                    { If AllowProtein=True then do not filter out the protein bases. Else keep only DNA (and ambiguity) bases }
    ExtraSeparator   : Boolean;                                                                    { Separate the sequences with an additional empty line }
    NoComments       : Boolean;                                                                    { Do not add comments AT ALL when saving file to disk. Just put one sequence on each line. Useful to build compact FASTA files. }

    constructor Create(aLog: TRamLog);
    destructor Destroy; override;
    procedure  Clear;
    procedure  FreeAllObjects;

    procedure  Add(Fasta: TFastaObj);                                          overload;
    procedure  Add(FullFileName: string; Bases: BaseString; Comments: string); overload;
    procedure  SetCapacity(Total: Integer);                                                        { Prealocate space }

    procedure  Save;
    procedure  SaveToFile(CONST FullFileName: string);
    function   LoadFromFile (CONST FullFileName: string; CleanEnds: Boolean= TRUE): Boolean;
    procedure  Split2Fasta;                                                                        { Split this multiFasta file in multiple individual FASTA files. OutputFolder= the path where the files were saved. }

    function   Count: Integer;
    property   Fasta[Index: Integer]: TFastaObj read getFasta;

    function   GetAllSequences: BaseString;
    procedure  UseNameAsComment;                                                                   { Replace original comments with the name of the FASTA file }
 end;


 
IMPLEMENTATION
USES ccIO;




{--------------------------------------------------------------------------------------------------
                              TMultiFasta
--------------------------------------------------------------------------------------------------}
constructor TMultiFasta.Create(aLog: TRamLog);
begin
 inherited Create;
 RamLog:= aLog;
 List:= TObjectList.Create;
 List.OwnsObjects:= TRUE;                                                                          { IT WILL FREE THE OBJECTS ON DESTROY! }
 
 Clear;

 RowLength        := ctDefaultRowLengt;
 ExtraSeparator   := FALSE;                                                                        { Separate the sequences with an additional empty line }
 NoComments       := FALSE;                                                                        { Do not add comments AT ALL when saving file to disk. Just put one sequence on each line. Useful to build compact FASTA files. }
 AllowProtein     := FALSE;                                                                        { If AllowProtein=True then do not filter out the protein bases. Else keep only DNA (and ambiguity) bases }
 Name_UseOrigName := TRUE;
 Name_UseIncrement:= FALSE;
 Name_UseComments := FALSE;
 Name_CustomName  := '';

 CleanedVectors   := '';
 OutputFolder     := '';                                                                           { The folder where the output files will be saved when the Split2Fasta function is used }
end;


destructor TMultiFasta.Destroy;
begin
 FreeAndNil(List);                                                                                 { Free also all contained FASTA objects }
 inherited Destroy;
end;


procedure TMultiFasta.Clear;
begin
 FileName:= '';
 List.Clear;                                                                                       { Asta elibereaza obiectele pt ca OwnsObjects=TRUE }
end;


procedure TMultiFasta.FreeAllObjects;                                                              { Remember List.OwnsObjects is by default TRUE so I don't have to call FreeAllObjects manually. }
VAR i: Integer;
    Fasta: TObject;
begin
  for i:= List.Count-1 downto 0 DO
   begin
    Fasta:= List[i];                                                                               { Get object }
    List.Extract(Fasta);                                                                           { Remove it from list }
    FreeAndNil(Fasta);                                                                             { Free it }
   end;
end;


procedure TMultiFasta.SetCapacity(Total: Integer);                                                 { Prealocate space }
begin
 List.Capacity:= Total
end;



procedure TMultiFasta.Add(Fasta: TFastaObj);
begin
 List.Add(Fasta);
end;


procedure TMultiFasta.Add(FullFileName: string; Bases: BaseString; Comments: string);
VAR aFasta: TFastaObj;
begin
 aFasta:= TFastaObj.Create(RamLog);                                                                { Nu ii dau aici FREE. De asta se va ocupa 'List' }
 aFasta.FileName:= FullFileName;
 aFasta.RowLength:= RowLength;                                                                     { Default row length for FASTA file. After this length the BASES string will be broken to a new row }
 aFasta.Comment  := Comments;
 aFasta.Bases    := Bases;
 Add(aFasta);
end;







{--------------------------------------------------------------------------------------------------
   LOAD/SAVE
--------------------------------------------------------------------------------------------------}
function TMultiFasta.LoadFromFile(CONST FullFileName: string; CleanEnds: Boolean= TRUE): Boolean;       { If CleanEnd = True, then clean up sequences during import: Remove the N bases at the ends of the sequence (not also at the middle) }
VAR Lines: TStringList;
    i, CurLine, TotalComments: integer;
    Baze: BaseString;
    sComment: string;
    aFasta: TFastaObj;
begin
 Clear;                                                                                                 { Golesc variabilele de continutul anterior }
 Result:= FALSE;
 CurLine:= 0;

 if NOT FileExists(FullFileName) then
  begin
   RamLog.AddError('The file does not exist'+CRLF+ FullFileName);
   EXIT(FALSE);
  end;

 if NOT IsPlainText(FullFileName) then
  begin
   RamLog.AddError('Only the following Fasta formats are supported: fasta, fa, fas, fst, fsa, fsta, seq, txt.'+CRLF+ FullFileName);
   EXIT(FALSE);
  end;

 FileName:= FullFileName;
 TRY

   Lines:= TStringList.Create;
   TRY
    { Open in SHARE mode }                                                                       { TRICK! Deschid fisierul cu StringFromFileA care foloseste 'fmShareDenyNone' ca sa nu am share violation cand fisierul e deschis si in alta aplicatie, cum ar fi Word }
    Lines.Text:= String(StringFromFileA(FullFileName));                                          { Am avut un EFOpenError report aici. Cannot open file "C:\Users\ahartmann\Documents\HCP2009\plasmide\X3 allcontigs.txt". Probably locked. }

    { Delete empty lines }
    for i:= Lines.Count-1 downto 0 DO
      if Lines[i]= ''
      then Lines.Delete(i);

    { Empty file? }
    if Lines.Count = 0 then
     begin
      RamLog.AddError('This FASTA file is blank');
      EXIT(FALSE);
     end;

    { Make it work with Seq/Txt files, that have no comment }
    if Pos('>', Lines.text) = 0
    then Lines.Insert(0, '>');

    { Protection against FASTA files that have no comment line }
    for i:= 0 to Lines.Count-1 DO
      if Pos('>', Lines.Strings[i]) > 0 then                                                { Search the first ROW containig a valid entry (comment) }
       begin
        CurLine:= i;
        Break;
       end;

    TotalComments:= CountAppearance('>', Lines.Text);
    if TotalComments< 1
    then RamLog.Addhint('No comment line detected!'); { Seq and txt files have no comment! }

    { Split MULTI-FASTA in its sub-components }
    REPEAT

     { Detect new comment }
     if Pos('>', Lines.Strings[CurLine]) < 1 then                                         { Daca nu am dat de semnul '>' inseamnca ca blockurile FASTA sunt sparate de o linie goala }
      begin
        inc(CurLine);                                                                     { Treci la urmatoarea linie }
        if CurLine >= Lines.Count
        then EXIT                                                                           { End of file. Exit. }
        else Continue;                                                                      { Jump to next line, maybe it contains a '>' sign }
      end;

     sComment:= Lines.Strings[CurLine];

     { After getting the comments, go to the next line and pick up the bases }
     inc(CurLine);

     { Build bases }                                                                        { Current sample may have bases spread over one or more lines. Search all rows that contains those bases, until I find the next comment line }
     Baze:= '';
     WHILE (CurLine< Lines.Count)                                                         { pana nu mai am linii }
      AND (Pos('>', Lines.Strings[CurLine])< 1) DO                                        { pana dau de un comentariu }
       begin
        Baze:= Baze+ Lines.Strings[CurLine];
        inc(CurLine);
       end;

     { Has ads? }
     if AdsDetected(Baze)
     then AllowProtein:= TRUE;                                                              { In Demo versions I will add a nag text. If the filter is active, it will remove that text so I have to disable the filter. }

     { FILTER invalid bases }
     if CleanEnds
     then Baze:= TrimSequenceNs(Baze);                                                      { Remove the N bases from the ends of a sequences }

     { Create new FASTA object }
     aFasta:= TFastaObj.Create(RamLog);                                                     { Nu ii dau aici FREE. De asta se va ocupa 'List' }
     TRY
       aFasta.FileName    := FullFileName;                                                  { See: TSmplImporter.GenerateVirtualNames }
       aFasta.AllowProtein:= AllowProtein;
       aFasta.RowLength   := RowLength;                                                     { Default row length for FASTA file. After this length the BASES string will be broken to a new row }
       aFasta.Comment     := sComment;
       aFasta.Bases       := Baze;
       aFasta.IsPart      := TotalComments> 1;                                              { Shows if this GBK is part of a multi-GBK object }
       aFasta.ParentType  := bfFAS;

       { Store it }
       if aFasta.NoOfBases >= ctMinimimSequence                                                              { If I have enough bases }
       then List.Add(aFasta)
       else
          begin
           RamLog.AddError('A sub-sequence was not loaded because it is too short!');
           FreeAndNil(aFasta);
          end;

     EXCEPT
       RamLog.AddError('Sub-sequences is invalid.');
       FreeAndNil(aFasta);
     END;
   UNTIL CurLine >= Lines.Count-1;

   { Virtual names }
   if TotalComments> 1
   then GenerateVirtualNames;

   Result:= List.Count> 0;
  FINALLY
    FreeAndNil(Lines);
  END;

 EXCEPT
    on EOutOfResources DO
     begin
      Result:= FALSE;
      RamLog.AddError('Not enough resources available to process this MultiFasta file.'+ CRLF+ FullFileName);
     end;

      on E: Exception DO
       begin
        Result:= FALSE;
        RamLog.AddError('Cannot load sample. '+ E.Message);                                            { Excpetion is handled here and does not propagate anymore }
       end;
 END;
end;



procedure TMultiFasta.SaveToFile(CONST FullFileName: string);
begin
 FileName:= FullFileName;
 Save;
end;



procedure TMultiFasta.Save;
begin
 StringToFileA(FileName, GetAllSequences, woOverwrite);                                            { StringToFile_A nu face "Test write access" }
end;



procedure TMultiFasta.UseNameAsComment;                                                            { Replace original comments with the name of the FASTA file }
VAR i: Integer;
begin
 for i:= 0 to Count-1 DO
  getFasta(i).Comment:= ExtractOnlyName(getFasta(i).FileName);
end;



procedure TMultiFasta.GenerateVirtualNames;                                                        { When I load a MultiFasta file, all cubs will have the same name, so here I generate a unique name for each cub }
VAR i: Integer;
    LeadingZeros: string;
begin
 for i:= 0 to Count-1 DO
  begin
   LeadingZeros:= '['+ LeadingZerosAuto( i2s(i), Count-1 )+ '] ';
   Fasta[i].FileName:= AppendBeforeName(Fasta[i].FileName, LeadingZeros);              { This will make from 'multi.fasta' something like '[01] multi.fsta', '[02] multi.fsta', '[03] multi.fsta' etc }
  end;          
end;











procedure TMultiFasta.Split2Fasta;                                                         { Split this multiFasta file in multiple individual FASTA files. OutputFolder= the path where the files were saved. }
CONST
   ctFastaExt= '.Fasta';
VAR i, iLeftChars: integer;
    sNameCounter, s, OutName, sNameOrig, sNameComment: string;
    aFasta: TFastaObj;
begin
 { Create path }
 if NOT ForceDirectoriesMsg(OutputFolder) then
  begin
   MesajError('Cannot create folder:'+ CRLF+ OutputFolder);
   Exit;
  end;
  
 Name_CustomName:= CorrectFilename(Name_CustomName, '_');

 { Enumerate all internal FASTA objects } 
 for i:= 0 to Count-1 do
  begin
   aFasta:= TFastaObj(List[i]);

   if Name_UseOrigName
   then sNameOrig:= ExtractOnlyName(aFasta.FileName);

   if Name_UseIncrement AND (Count> 1)                                                             { Applies only if the input file contains more than one subsamples }
   then sNameCounter:= ' '+ LeadingZerosAuto(i2s(i+1), Count);

   if Name_UseComments  then
    begin
     iLeftChars:= MAX_PATH - Length(OutputFolder+ sNameOrig+ Name_CustomName+ sNameCounter+ ctFastaExt)- 1;  { don't allow more than MAX_PATH (260) chars in path. Nu merge fara acel 1 de la coada }
     s:= Trim(aFasta.Comment);                                                                     { asta e obligatoriu zice Cristina pt ca vrea sa obtina nume de fiisier fara spatii la urma }
     s:= StringReplace(s, '>', '', [rfReplaceAll, rfIgnoreCase]);                                  { asta vine inainte de CorrectPath }
     s:= CorrectFilename(s, '_');
     s:= Trim(s);                                  
     sNameComment:= ' '+ system.COPY(s, 1, iLeftChars);
    end;   

   OutName:= OutputFolder+ Trim(sNameOrig+ Name_CustomName+ sNameComment+ sNameCounter+ ctFastaExt);
   aFasta.SaveAsFasta(OutName, FALSE);                                                          { This also tests for write access }
  end;
end;


function TMultiFasta.Count: Integer;
begin
 Result:= List.Count;
end;



function TMultiFasta.getFasta(Index: Integer): TFastaObj;
begin
 Result:= TFastaObj(List[Index]);
end;


function TMultiFasta.GetAllSequences: BaseString;
VAR i: Integer;
    aFasta: TFastaObj;
    Separator: String;
begin
 Result:= '';

 if ExtraSeparator                                                                                 { Separate the sequences with an additional empty line }
 then Separator:= LBRK
 else Separator:= CRLF;

 { Enumerate all internal FASTA objects }
 for i:= 0 to Count-1 do
  begin
    aFasta:= Fasta[i];
    if NoComments
    then Result:= Result+ aFasta.Bases+ Separator                                          { Do not add comments AT ALL when saving file to disk. Just put one sequence on each line. Useful to build compact FASTA files. }
    else Result:= Result+ aFasta.Comment+ CRLF+ aFasta.Bases+ Separator;
  end;
end;





END.
