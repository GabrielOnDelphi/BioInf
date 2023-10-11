# BioInf

Library for reading microbiology DNA files (ABI, SCF, FASTA, etc).
Requires LightSaber base library.

Has support for:
 * SNP (Single Nucleotide Polymorphism)
 * QV (base Quality Value)
 * Integrated trimming engine (automatic bad-end trimming)
 * Recalling of bad (N) peaks using proprietary algorithm 
 * Automatic sequence direction (F/R)
 * Reading and displaying the chromatogram data
 * read/write SCF files  
 * read ABI/AB/AB1/AB! files  
 * read/write FASTA files  
 * read/write GeneBank (GBK) files.  
  
The import functionality is accessible with a call to a single function: TCubeImport.Import(FileName)  
  
  TCube = class(TCubeAbstractSnp)    
  public   
   {IMPORT}   
   function  Import(CONST FileName: string): Boolean;                                                         { Fasta and GBK files are not supported because they may contain more than one sample! }  
   function  AssignGBK     (CONST Gbk : TGbkObj): Boolean;  
   function  AssignScf     (CONST SCF : TScfObj): boolean;  
   function  AssignAbi     (CONST ABI : TAbiObj; EditedField: boolean= True): Boolean;                      
   {EXPORT}                           
   function  ExportAs      (CONST aFileName: string; Clean: Boolean; OUT BasesLeft: Integer): Boolean;          
   procedure ExportAsFASTA (CONST FileName : string; Clean: Boolean; OUT BasesLeft: Integer); overload;   
   procedure ExportToSCF   (CONST aFileName: string; Clean: Boolean);            overload;                    { This will remove all GAPs }  
   function  ExportToSCF   (CONST                    Clean: Boolean): string;    overload;                    { Returns the name used to save the file }  
   function  AsFasta       (CONST aFileName: string; Clean: Boolean): TFastaObj; overload;                    { Return a TFastaObj object built from this cube }  
  
   procedure SaveChromaAsBmp(CONST Nume: string);  
 end; 

 The code was tested on millions of DNA samples over the years and works flawlessly. 

-----------

Next to come:   
Code to import NextGen DNA sequence files such as SFF. SFQ, FAS.

  Features:
      + TrimEngine
      + Assign(TSample)
      + Average Quality
      + Highest QV
      + FindAdaptor
      + Clip Left/Right
      + Good Bases
      + can work with sequences with length over 65535 bases
      + Implements a single SFF read
      + NoOfGoodBases      
      + GoodBases
      + Conversion to Fasta
      + BuildQVGraph   
      etc
      
