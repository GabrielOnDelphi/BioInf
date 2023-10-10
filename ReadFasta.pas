
UNIT ReadFasta;

{==============================================================================================
   Heracle BioSoft SRL
   2016.07.13

   Single Fasta file
==============================================================================================}

INTERFACE
USES
   System.Classes, System.SysUtils, ccCore, ccINIFile, CubicDNA, clRamLog, ccRichLog, ReadBasicSample;

TYPE
  TFastaObj = class(TBasicSample)
   private
     procedure setBases (CONST Value: BaseString);
   protected
   public
     Vectors     : RVectorData;
     RowLength   : integer;
     AllowProtein: Boolean;                                                             { If AllowProtein=True then do not filter out the protein bases. Else keep only DNA (and ambiguity) bases }
     constructor Create(aLog: TRamLog);
     function   LoadFromFile (FullFileName: string; CleanEnds: Boolean= true): Boolean;
     procedure  SaveAsFasta  (CONST sFileName: string; Wrap: Boolean);                  { SAVE AS FASTA }
     procedure  SaveAsSEQ    (CONST sFileName: string);
     procedure  SaveAsTxt    (CONST sFileName: string);                                 { SAVE AS TXT }
     procedure  Save(Wrap: Boolean);
     procedure  Save_AutoDetect;                                                         { returns '' is everything was OK }

     property   Bases: BaseString read FBases     Write SetBases;
  end;



IMPLEMENTATION
USES ccIO, cmWrapString;








{--------------------------------------------------------------------------------------------------
   CONSTRUCTOR
--------------------------------------------------------------------------------------------------}

constructor TFastaObj.Create(aLog: TRamLog);                                                            { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
begin
 inherited Create(alog);                                                                                { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }
 AllowProtein:= FALSE;                                                                                  { If AllowProtein=True then do not filter out the protein bases. Else keep only DNA (and ambiguity) bases }
 RowLength := ctDefaultRowLengt;
end;




{--------------------------------------------------------------------------------------------------
   LOAD FROM FILE
--------------------------------------------------------------------------------------------------}
function TFastaObj.LoadFromFile(FullFileName: string; CleanEnds: Boolean= true): Boolean;               { If CleanEnd = True, then clean up sequences during import: Remove the N bases at the ends of the sequence (not also at the middle) }
VAR Lines: TStringList;
    Baze: string;
    i: Integer;
begin
 Clear;                                                                                                 { ABSOLUTELLY NECESAR }

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

 ParentType:= bfFAS;
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

    { Multi fasta? }
    IsPart:= CountAppearance('>', Lines.Text) > 1;
    if IsPart
    then RamLog.AddWarn('This is a multiFASTA file. Only the first sequence was loaded!');

    { Process comments }
    if Pos('>', Lines.Strings[0]) > 0
    then
       begin
        Comment:= Lines.Strings[0];
        Lines.Delete(0);
       end
    else
       GenerateComment;

    { Delete the rest of the sequences }
    if IsPart then
     begin
      i:= Pos('>', Lines.Text);                                                                  { Find the sencond comment }
      Lines.Text:= system.COPY(Lines.Text, 1, i-1);                                              { Copy from start to second comment }
     end;

    { Read bases }
    Baze:= Lines.Text;

    { Has ads? }
    if AdsDetected(Baze)
    then AllowProtein:= TRUE;                                                                    { In Demo versions I will add a nag text. If the filter is active, it will remove that text so I have to disable the filter. }

    { FILTER invalid bases }
    if CleanEnds
    then Baze:= TrimSequenceNs(Baze);                                                            { Remove the N bases from the ends of a sequences }   { Cazul Lopez. Clean the numbers and invalid characters from this string. }

    { Has bases? }
    if Length(Baze) < ctMinimimSequence then
     begin
      RamLog.AddError('The sequence is too short!');
      EXIT(FALSE);
     end;

    Bases:= Baze;
    Result:= TRUE;
   FINALLY
     FreeAndNil(Lines);
   END;

 EXCEPT
   on E: Exception DO
    begin
     Result:= FALSE;
     RamLog.AddError('Cannot load sample: '+ FullFileName+ ' '+ E.Message);                      { Excpetion is handled here and does not propagate anymore }
    end;
 END;
end;





{--------------------------------------------------------------------------------------------------
   SAVE TO FILE
--------------------------------------------------------------------------------------------------}
procedure TFastaObj.SaveAsFasta(CONST sFileName: string; Wrap: Boolean);                           { SAVE AS FASTA }
VAR sComment, sBody, sBases: String;
begin
 FileName:= sFileName;

 { COMMENTS }
 if CommentIsEmpty                                                                                 { daca nu am nici un comentariu atunci creez eu unul }
 then GenerateComment;

 if Wrap
 then sBases:= cmWrapString.WrapStringForced(Bases, RowLength)
 else sBases:= Bases;
 sComment:= Comment;                                                                               { RULE: The TFastaObj.SaveToDisk is responsible for mixing the Comments and MetaData in one single string }
 sBody:= sComment+ CRLF+ sBases;
 StringToFileA(FileName, sBody, woOverwrite);                                                      { StrinToFileA face si "Test write access" }
end;



procedure TFastaObj.SaveAsSeq(CONST sFileName: string);                                            { SAVE AS SEQ }
begin
 ParentType:= bfSEQ;                                                                               { < aici nu trebuie sa trunchez liniile la 80 caractere? }
 StringToFileA(sFileName, BASES, woOverwrite);                                                     { StrinToFileA face si "Test write access" }
end;



procedure TFastaObj.SaveAsTxt(CONST sFileName: string);                                            { SAVE AS TXT }
begin
 ParentType:= bfTXT;
 StringToFileA(sFileName, BASES, woOverwrite);                                                     { StrinToFileA nu face  "Test write access" }
end;



procedure TFastaObj.Save_AutoDetect;                                                               { Auto detect save type }
begin
 if IsFas(FileName)
 then SaveAsFasta(FileName, TRUE)
 else
    if IsSeq(FileName)
    then SaveAsSeq(FileName)
    else
      if IsTXT(FileName)
      then SaveAsTxt(FileName)
      else RamLog.AddError('Cannot autodetect file type: '+ FileName);
end;



procedure TFastaObj.Save(Wrap: Boolean);                                                           {TODO: Later. Rename it back to 'Save' }
begin
 SaveAsFasta(FileName, Wrap);
end;


procedure TFastaObj.setBases(CONST Value: BaseString);
begin
 FBases:= CubicDNA.CleanSequence(Value, FALSE);
 FBases:= CubicDNA.RemoveGaps(FBases);
end;



END.
