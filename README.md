# BioInf

Library for reading microbiology DNA files (ABI, SCF, FASTA, etc).
Requires LightSaber base library.


Allows you to:  
 * read/write SCF files  
 * read ABI/AB/AB1/AB! files  
 * read/write FASTA files  
 * read/write GeneBank (GBK) files.  
  
 All this functionality is accessible with simple access to TCubeImport.Import(FileName)  
  
  TCube = class(TCubeAbstractSnp)    
  public   
   {IMPORT}   
   function  Import(CONST FileName: string): Boolean;                                                         { Fasta and GBK files are not supported because they may contain more than one sample! }  
   function  AssignGBK     (CONST Gbk : TGbkObj): Boolean;  
   function  AssignScf     (CONST SCF : TScfObj): boolean;  
   function  AssignAbi     (CONST ABI : TAbiObj; EditedField: boolean= True): Boolean;                      
   function  AssignSecv    (CONST SECV: TFastaObj): Boolean;  
   property  Version: string  read FParentVer  Write FParentVer;                                            
  
   {EXPORT}  
   function  ExportBack    (CONST Clean: Boolean; OUT BasesLeft: Integer): Boolean;                            
   function  ExportAs      (CONST aFileName: string; Clean: Boolean; OUT BasesLeft: Integer): Boolean;         
   procedure ExportAsFASTA (CONST                    Clean: Boolean; OUT BasesLeft: Integer); overload;        
   procedure ExportAsFASTA (CONST FileName : string; Clean: Boolean; OUT BasesLeft: Integer); overload;   
   procedure ExportAsSecv  (CONST NewName  : string; Clean: Boolean; OUT BasesLeft: Integer);  
   procedure ExportToSCF   (CONST aFileName: string; Clean: Boolean);            overload;                    { This will remove all GAPs }  
   function  ExportToSCF   (CONST                    Clean: Boolean): string;    overload;                    { Returns the name used to save the file }  
   function  AsFasta       (CONST aFileName: string; Clean: Boolean): TFastaObj; overload;                    { Return a TFastaObj object built from this cube }  
   function  AsFasta       (CONST                    Clean: Boolean): TFastaObj; overload;                    { Same as above but the user doesnt have to provide a file name. It is automatically assumed from Cube's name }  
  
   procedure SaveChromaAsBmp(CONST Nume: string);  
 end; 

 The code was tested on millions of DNA samples over the years and works flawlessly. 
