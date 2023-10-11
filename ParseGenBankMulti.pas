UNIT ParseGenBankMulti;

{=============================================================================================================
 Gabriel Moraru
 2016.07
==============================================================================================================

 1. Load the whole text in FRawText.
 2. Find how many GBK files I have incuded in this file
 3. Split the raw text in its sub-components (GBK files)
 4. Add each sub-GBK to the GbkList of objects (GbkList)
 5. Parse each object.

===============================================================================}

INTERFACE
USES
   System.SysUtils, System.Classes, System.Contnrs,
   ccCore, ccINIFile, clRamLog, ccRichLog, CubicDNA, ReadBasicSample, ParseGenBankR, ParseGenBankW;

TYPE
 TMultiGbk = class(TObject)                                                                             { TObject is never directly instantiated. Although it does not use programming language features that prevent instantiation, TObject is an abstract class. }
   private
     FRawText  : TStringList;
     FFileName : string;
     function  getBank (Index: Integer): TGbkObjExp;
   protected
     RamLog: TRamLog;
     GbkList: TObjectList;                                                                              { List of GBK objects }
     procedure generateVirtualNames;                                                                    { When I load a MultiFasta file, all cubs will have the same name, so here I generate a unique name for each cub }
   public
     constructor Create(aLog: TRamLog);
     destructor Destroy; override;
     procedure  Clear;
     {}
     function LoadFromFile (FullFileName: string): Boolean;
     function Count: Integer;
     property FileName           : string      Read FFileName  Write FFileName;
     property GBK[Index: integer]: TGbkObjExp  Read getBank;  default;
 end;



IMPLEMENTATION

USES
   ccIO;



constructor TMultiGbk.Create;
begin
 inherited Create;                                                                                      { Should I call "inherited" in the constructor of a class derived from TObject or TPersistent? Yes. It does nothing, true, but it's harmless. I think there is value in being consistent about always calling the inherited constructor, without checking to see if there is, in fact, an implementation. Some will say that it's worth calling inherited Create because Embarcadero might add an implementation for TObject.Create in the future, but I doubt this is true; it would break existing code which does not call inherited Create. Still, I think it is a good idea to call it for the reason of consistency alNone. }

 RamLog:= aLog;
 FRawText:= TStringList.Create;
 GbkList:= TObjectList.Create;
 GbkList.OwnsObjects:= TRUE;                                                                            { IT WILL FREE THE OBJECTS ON DESTROY! }

 Clear;                                                                                                 { Asta e intotdeauna ultimul }
end;



destructor TMultiGbk.Destroy;
begin
 FreeAndNil(FRawText);
 FreeAndNil(GbkList);

 inherited;
end;



procedure TMultiGbk.Clear;
begin
 FFileName:= '';
 FRawText.Clear;
 GbkList.Clear;                                                                                         { Asta elibereaza obiectele pt ca OwnsObjects=TRUE }
end;




function TMultiGbk.Count: Integer;
begin
 Result:= GbkList.Count;
end;




procedure TMultiGbk.generateVirtualNames;                                                               { When I load a MultiGbk file, all cubs will have the same name, so here I generate a unique name for each cube }
VAR i: Integer;
    LeadingZeros: string;
begin
 if GbkList.Count > 1 then
  for i:= 0 to Count-1 DO
   begin
    LeadingZeros:= '['+ LeadingZerosAuto( i2s(i), Count-1 )+ '] ';
    getBank(i).FileName:= AppendBeforeName(getBank(i).FileName, LeadingZeros);                          { This will make from 'multi.fasta' something like '[01] multi.fsta', '[02] multi.fsta', '[03] multi.fsta' etc }
   end;
end;




function TMultiGbk.getBank(Index: integer): TGbkObjExp;
begin
 Assert(Index< GbkList.Count);
 Result:= TGbkObjExp(GbkList[Index]);
end;










{--------------------------------------------------------------------------------------------------
   LOAD
--------------------------------------------------------------------------------------------------}
function TMultiGbk.LoadFromFile(FullFileName: string): Boolean;
VAR
   GbkStart, GbkEnd, CopyLine, i: Integer;
   GBK: TGbkObjExp;
begin
 RamLog.AddVerb('Loading: '+ FullFileName);

 { CHECK NAME }
 if NOT FileExists(FullFileName) then
   begin
     RamLog.AddError('GBK file does not exist! '+ FullFileName);
     EXIT(FALSE);
   end;

 if NOT IsGBK(FullFileName) then
   begin
     RamLog.AddError('This is not a GBK file: '+ FullFileName);
     EXIT(FALSE);
   end;

 { INIT }
 Clear;                                                                                        { ABSOLUT NECESAR }
 GbkEnd:= -1;
 FFileName:= FullFileName;

 TRY
   { LOAD }
   FRawText.LoadFromFile(FullFileName);

   { Break Multi-GBK into smaller pieces }
   GbkStart:= 0;
   for i:= 0 TO FRawText.Count-1 DO                                                            { search through all lines }
    if Pos(ctGbkEndOfSample, FRawText.Strings[i])= 1 then                                      { search separator -> AFLU UNDE TERMINA GENBANK-ul CURENT }
     BEGIN
      GbkEnd:= i;

      GBK:= TGbkObjExp.Create(RamLog);                                                         { FREED BY: 'GbkList.OwnsObjects:= TRUE' or below if the file is invalid }
      GBK.FileName:= FullFileName;

      { Add raw data }
      for CopyLine:= GbkStart to GbkEnd                                                        { includ si separatorul }
       DO GBK.addRawData(FRawText[CopyLine]);

      { Process raw data and load info into the appropriate fields }
      GBK.Parse;

      { Store current GBK }
      GbkList.Add(GBK);

      GbkStart:= GbkEnd+1;
     END;
 EXCEPT
   RamLog.AddError('GBK file cannot be loaded: '+ FullFileName);                               { Don't crash during import! }
 END;

 if GbkEnd= 0
 then RamLog.AddWarn('Cannot find end of the GBK file! '+ FullFileName);

 for i:= 0 to Count-1
   DO getBank(i).IsPart:= GbkList.Count > 1;                                                   { If multiGBK then mark all objects as so }

 if GbkList.Count = 0
 then RamLog.AddError('GBK file is invalid or empty! '+ FullFileName)
 else generateVirtualNames;                                                                    { When I load a MultiFasta file, all cubes will have the same name, so here I generate a unique name for each cub }

 Result:= GbkList.Count > 0;
end;



end.
