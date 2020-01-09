{
 Suplementary Information for

 Drivers of Geographic patterns of North American Language Diversity

 Marco Túlio Pacheco Coelho, Elisa Barreto Pereira,
 Hannah Haynie, Thiago F. Rangel, Patrick Kavanagh,
 Kathryn R. Kirby, Simon J. Greenhill, Claire Bowern,
 Russell D. Gray, Robert K. Colwell, Nicholas Evans,
 Michael C. Gavin

 Marco Túlio Pacheco Coelho
 Email: marcotpcoelho@gmail.com
 Michael C. Gavin
 Email: Michael.Gavin@colostate.edu

 This source code represents the carrying capacity simulation model
 published by Gavin et al. (2017) - Process‐based modelling shows
 how climate and demography shape language diversity.
 https://doi.org/10.1111/geb.1256

}




unit LibLPMMain;

interface

uses
  {$IfDef DLLShapes}LibFastShareMem,{$EndIf}
  SysUtils, Classes, Math, Linux, SyncObjs,
  {$IfNDef CONSOLE}
  ExtCtrls,
  {$EndIf}
  LibTypes, LibGIS, LibFiles, LibStats, LibPhyStats,
  LibGeometry, LibMCMC,
  LibMatrix, LibProbDistrib;

Type
   TUpdateProc = Procedure of object;

Type
   TModelFunction = (TMFPower, TMFLogistic, TMFExponential, EstimatedDensity);

Type
 TModelParameters = Record
   OutPutFolder, SimName: String;

   DataPath: String;
   MapShapeFile: String;
   MapShapeSmallResFile: String;
   EnvDataFile: String;
   SubsitModeParamFile: String;
   AreaParamFile: String;
   MaxPopSizeFile: String;
   NumInitLang: Integer;
   InitPopSize: Integer;
   PopPerLanguage: Integer;
   InitCell: Integer;

   ModelFunction: TModelFunction;

   CarryingCapParams: TDblVector;

   r: Double;                          // Intrinsic population growth rate
   NeighborDistance: Double;
   SpeciationTime: Integer;

   RandomSeed: Integer;

   ObsRichMap: TDblVector;
   ObsRSFD: TDblVector;
   ObsAvgRSMap: TDblVector;
  End;

Type
 TParComboRepsResults = Record
   SumRichPerCell: TDblVector;
   Sum2RichPerCell: TDblVector;
   SumAvRangeSizePerCell: TDblVector;
   Sum2AvRangeSizePerCell: TDblVector;
   AcumRSFD: TDblVector;
   NumLangs: TDblVector;
   RepCount: Integer;
  End;

 TModelReScaleResults = Record
   AvRangeSizePerCell: TDblVector;
   RichPerCell: TDblVector;
   RSFD: TDblVector;
  End;

 TModelFitStatsResults = Record
   nLangs: Double;
   RegResAvgRich: TRegResults;
   RegResAvgRS: TRegResults;
   KS: Double;
  End;


Type
 PLanguage = ^TLanguage;
 TPLangVec = Array of PLanguage;

 PMapCell = ^TMapCell;
 TPMapCellVec = Array of PMapCell;

 PLangMapCell = ^TLangMapCell;
 TPLangMapCellVec = Array of PLangMapCell;

 TColor = -$7FFFFFFF-1..$7FFFFFFF;

 TMapCell = Record
   IdCell: Integer;
   Coord: TPoint2D;

   Environment: TDblVector;
   Area: Double;
   Neighbors: TPMapCellVec;
   K: Integer;

   Languages: TPLangVec;
   TotPopSize: Integer;
   DeltaPop: Integer;
   PopSaturated: Boolean;

   ColonizedOnce: Boolean;
  end;

 TLangMapCell = Record
   IdCell: Integer;
   PosVec: Integer;

   MapCell: PMapCell;
   Language: PLanguage;

   K: Integer;
   LangPopSizeInCell: Integer;
   DeltaLangPop: Integer;

   PatchTag: Integer;
   PrevPatchTag: Integer;

   BorderCell: Boolean;
  End;
 TLangMapCellVec = Array of TLangMapCell;

 TLangPatch = Record
   PatchSeparationTime: TIntVector;   // Temporal distance to other patches
   Cells: TPMapCellVec;               // Cells that belong to the patch
  End;

 TLanguage = Record
   IdLang: Integer;

   OccupCells: TLangMapCellVec;
   SubsistModeProp: TDblVector;
   TotPopSize, MaxPopSize: Integer;
   DeltaPop: Integer;

   ExtantLang: Boolean;
   NeverExisted: Boolean;
   ExtinctLang: Boolean;
   RecentLang: Boolean;
   TransformedLang: Boolean;

   DaughterOf: Integer;
   AncestorOf: TPLangVec;
   BornAt: Integer;
   DiedAt: Integer;

   Perimeter: Cardinal;     // nCells_Border
   Area: Cardinal;          // nCells_Total

   NumPatches: Integer;
   Patches: Array of TLangPatch;

   ChangedThisStep: Boolean;
   ChangedLastStep: Boolean;
   PopSizeAccounted: Boolean;

   Color: TColor;
  End;

 TMaxPopClass = Record
   Min, Max: Double;          // Min and Max limits of pop size in each class

   ObsCount: Integer;          // Number of empirical languages that fall in the class
   ObsFrequency: Double;      // Frequency of languages that fall in the class

   SimCount: Integer;          // Number of simulated languages that fall in the class
   SimFrequency: Double;      // Frequency of simulated languages that fall in the class

   MaxPopVec: TDblVector;     // Vector with observed (empirical) language pop sizes
  End;

 TMaxPopDistrib = Array of TMaxPopClass;

 TSim = Class(TThread)
   Private
     SimId: Integer;

     // Random number generator and sample function is here to assure that results can be replicated
     LstRandSeed: Integer;           // Random Seed at the begining of the simulation
     RandSeed: Integer;              // Current Seed
     Procedure Randomize;            // Create a new random seed
     Function Random: Extended;      // Generates a random number
     Function Sample(List: TIntVector; n: Integer = -1; WithReplacement: Boolean = False): TIntVector;
   Public
    DebugMode: Boolean;

    Halted: Boolean;
    StopRepWhenHalt: Boolean;
    RunStepsContinuosly: Boolean;
    RunRepsContinuosly: Boolean;

    ModelParameters: TModelParameters;

    MapAnim: Boolean;
    GIFAnim: Boolean;
    AVIAnim: Boolean;
    MapAnimSampInt: Integer;
    MapAnimExpInt: Integer;
    GraphAnim: Boolean;
    GraphAnimSampInt: Integer;
    GraphAnimExpInt: Integer;
    AnimDelay: Double;
    PlotPhy: Boolean;

    PrjName: String;
    OutPutFolder: String;
    FileOutReplicate: ^TextFile;
    FileOutReplicateWrite: ^RTL_CRITICAL_SECTION;

    ModelSet: Boolean;
    AutoProg: Boolean;

    MapShape: TShapeFile;
    MapShapeSmallRes: TShapeFile;
    MapCells: Array of TMapCell;        // All cells in the map
    Langs: Array of TLanguage;          // All langauges

    LangSeq: TIntVector;
    CellSeq: TIntVector;

    WrkCells: TIntVector;
    WrkLangs: TIntVector;

    MaxPopVec: TDblMatrix;
    MaxPopDistrib: TMaxPopDistrib;

    InitCell: Integer;

    TotNumTimeSteps: Integer;
    TotNumReps: Integer;
    TotNumCells: Integer;
    TotNumLangs: Integer;
    TotNumLangsEver: Integer;
    TotNumEnvVars: Integer;
    TotNumCarryingCapParams: Integer;
    TotNumCellsColonizedOnce: Integer;

    ActNumReps: Integer;

    TotNumPatches: Integer;

    TimeStep: Integer;
    Rep: Integer;
    SimSpeed: Double;

    ErrMsg: String;

    ReScaledRep: TModelReScaleResults;
    FitStatsRep: TModelFitStatsResults;
    AvgFitStats: TModelFitStatsResults;
    FitIndex: Double;

    ParComboResults: TParComboRepsResults;

    RunNextStepEvent: TEvent;
    RunNextRepEvent: TEvent;

    StepFinished: TEvent;
    RepFinished: TEvent;
    SimFinished: TEvent;

    UpdateNextStepExternalProcedure: TUpdateProc;


    Constructor Create(CreateSuspended: Boolean;
                       pStopRepWhenHalt: Boolean = False;
                       pRunStepsContinuosly: Boolean = False;
                       pRunRepsContinuosly: Boolean = False);

    Procedure RunNextStep;

    Function SetUpModel(ReSetMapData: Boolean = True): Boolean;
    Function ReSetModel: Boolean;
    Procedure ClearParComboResults;

    Function CellChangePop(Lang: PLanguage; MapCell: PMapCell; DeltaNIndiv: Integer): PLangMapCell; overload;
    Function CellChangePop(LangMapCell: PLangMapCell; DeltaNIndiv: Integer): PLangMapCell; overload;

    Function CellChangeDeltaPop(Lang: PLanguage; MapCell: PMapCell; DeltaNIndiv: Integer): PLangMapCell; overload;
    Function CellChangeDeltaPop(LangMapCell: PLangMapCell; DeltaNIndiv: Integer): PLangMapCell; overload;

    Procedure RemoveFromWrkCells(MapCell: TMapCell);
    Procedure RemoveFromWrkLangs(Language: PLanguage);

    Procedure AddtoWrkCells(MapCell: TMapCell);
    Procedure AddToWrkLangs(Language: PLanguage);

    Function GetPopSize(Lang: PLanguage; MapCell: PMapCell): Integer;       overload;
    Function GetPopSize(LangMapCell: PLangMapCell): Integer;                overload;
    Function GetDeltaPopSize(Lang: PLanguage; MapCell: PMapCell): Integer;  overload;
    Function GetDeltaPopSize(LangMapCell: PLangMapCell): Integer;           overload;

    Function GetAllNeighbors(LangCell: PLangMapCell): TPMapCellVec;
    Function GetNeighborsDifferentLanguage(LangCell: PLangMapCell): TPLangMapCellVec;
    Function CountNeighborsSameLanguage(LangCell: PLangMapCell): Integer;
    Function GetNeighborsSameLanguage(LangCell: PLangMapCell): TPLangMapCellVec; overload;
    Function GetNeighborsSameLanguage(Language: PLanguage; Neighbors: TPMapCellVec): TPLangMapCellVec; overload;

    Function CreateNewLanguage(AncestorLang: PLanguage; SubsistModeProp: TDblVector): PLanguage; overload;
    Function CreateNewLanguage(AncestorLang: PLanguage; Cell: PMapCell; SubsistModeProp: TDblVector): PLanguage; overload;
    Function CreateNewLanguage(AncestorLang: PLanguage; CellOccup: TPMapCellVec; SubsistModeProp: TDblVector): PLanguage; overload;
    Function ChangeLanguageIdentity(AncestorLang: PLanguage): PLanguage;

    Function CalcCarryingCapacity(Lang: PLanguage; MapCell: PMapCell): Integer; overload;
    Function CalcCarryingCapacity(LangMapCell: PLangMapCell): Integer;          overload;
    Function CalcCarryingCapacity(MapCell: PMapCell): Integer;                  overload;

    Function SampleMaxPopSize: Integer;
    Procedure UpdatePopSizeFrequency(Language: PLanguage);

    Function GetLangMapCell(Lang: PLanguage; MapCell: PMapCell): PLangMapCell;
    Function GetLangMapCellVec(Lang: PLanguage; MapCellVec: TPMapCellVec): TPLangMapCellVec;
    Function IntRND(Min, Max: Double): Integer;

    Function ReScale: TModelReScaleResults;
    Function CalcRepFitStats: TModelFitStatsResults;
    Procedure CalcAvgRepFitStats;
    Procedure StoreRepResults;
    Function CalcRepFitIndex: Double;
    Function CalcThreadFitIndex: Double;

    Procedure ExportParComboEstResults;

    Function WritePhy(Parent: Integer; t: Integer): String;
    Procedure ExportResults;
  protected
    procedure Execute; override;
    Procedure UpdateNextStep;
    Procedure UpdateNextRep;
    Procedure UpdateSimFinished;
    Procedure ShowError;
  End;

 TSimParallel = Class
   ModelParameters: TModelParameters;
   GibbsSamp: TGibbs;

   SimVec: Array of TSim;
   HandleVec: Array of THandle;

   FileOutReplicate: TextFile;
   FileOutReplicateWrite: TRTLCriticalSection;
   FileOutSim: TextFile;

   AvgFitStats: TModelFitStatsResults;
   FitIndex: Double;
   SimFinished: TEvent;

   MapShape: TShapeFile;
   MapShapeSmallRes: TShapeFile;
   MapCells: Array of TMapCell;        // All cells in the map
   CellSeq: TIntVector;
   TotNumCells: Integer;
   TotNumEnvVars: Integer;

   MinLinParValSearch, MaxLinParValSearch: Double;
   MinQuadParValSearch, MaxQuadParValSearch: Double;

   BestFitLinParVal: Double;
   BestFitLinParFitIndex: Double;
   BestFitQuadParVal: Double;
   BestFitQuadParFitIndex: Double;

   private
    // Only used by optimization/search functions
    tnReps, tnThreads: Integer;

   public
     Function Execute(ModelParameters: TModelParameters;
                      nReps: Integer = 120;
                      nThreads: Integer = -1): Double;

     Function ExecToGibbsSamp(ProposalParams: TDblVector): Double;
     Procedure SearchGibbs(ModelParameters: TModelParameters;
                           DisturbFuncs: TDisturbFuncs;
                           var MCMCChain: TChain;
                           nReps: Integer = 120;
                           nThreads: Integer = -1);
  End;

Type
  TPatchCellItems = Record
    Found: Boolean;
    IdCell: Integer;
    AltPatchNum: Integer;
    AltPatchPos: Integer;
    Link: PMapCell;
   End;

  TOnePatchControl = Record
    PatchesLostCellsFrom: TIntVector; // List of patches that had cells now in this patch
    PatchCells: Array of TPatchCellItems;
   End;

  TPatchControl = Array of TOnePatchControl;

Type
  TLangPatchInfo = Record
    PatchSeparationTime: TIntVector;   // Temporal distance to other patches
    Cells: TPMapCellVec;
   end;


{$IfDef CONSOLE}
Procedure ErrorMessage(Msg: String);
{$EndIf}

implementation

{$IfNDef CONSOLE}
Uses
  UntControl, UntMap;
{$EndIf}

Constructor TSim.Create(CreateSuspended: Boolean;
                        pStopRepWhenHalt: Boolean = False;
                        pRunStepsContinuosly: Boolean = False;
                        pRunRepsContinuosly: Boolean = False);
 begin

  inherited Create(CreateSuspended);
  Self.FreeOnTerminate:= False;

  Halted:= False;

  MapAnim:= True;
  GIFAnim:= False;
  AVIAnim:= True;
  GraphAnim:= False;
  PlotPhy:= True;
  ModelSet:= False;
  AutoProg:= False;
  TimeStep:= 0;
  SimSpeed:= 0.001;

  TotNumReps:= 1;
  TotNumTimeSteps:= MaxInt;
  Rep:= 0;

  Self.StopRepWhenHalt:= pStopRepWhenHalt;
  Self.RunStepsContinuosly:= pRunStepsContinuosly;
  Self.RunRepsContinuosly:= pRunRepsContinuosly;

  RunNextStepEvent:= TEvent.Create;
  RunNextStepEvent.ResetEvent;

  RunNextRepEvent:= TEvent.Create;
  RunNextRepEvent.ResetEvent;

  StepFinished:= TEvent.Create;
  StepFinished.ResetEvent;

  RepFinished:= TEvent.Create;
  RepFinished.ResetEvent;

  SimFinished:= TEvent.Create;
  SimFinished.ResetEvent;

 end;

Procedure TSim.Execute;
 begin
  // If any replicate fails, the actual number of replicates is reduced
  ActNumReps:= TotNumReps;

  ClearParComboResults;

  // Simulation Loop
  While True do
   begin
    Rep:= Rep + 1;

    PrjName:= 'Rep' + IntToStr(Rep);

    If Rep > TotNumReps then
      Break;

    If Finished or Terminated then
      Break;

    ReSetModel;

    // Steps loop
    While True do
     begin
      // Waiting for the command to run the next step
      If (not RunStepsContinuosly) then
       begin
        RunNextStepEvent.ResetEvent;

        // Update status of the simulation here
        UpdateNextStep;

        WaitForSingleObject(RunNextStepEvent.Handle, INFINITE);
       end;

      If Finished or Terminated then
        Break;

      // Run the next step here.
      RunNextStep;

      StepFinished.SetEvent;
      StepFinished.ResetEvent;

      If Finished or Terminated then
        Break;

      // Continuous run until halt?
      If (StopRepWhenHalt and Halted) or (TimeStep > TotNumTimeSteps) then
        Break;
     end;

    RepFinished.SetEvent;
    RepFinished.ResetEvent;

    If Finished or Terminated then
      Break;

    ReScaledRep:= ReScale;
    FitStatsRep:= CalcRepFitStats;
    CalcRepFitIndex;
    CalcAvgRepFitStats;
    StoreRepResults;

    UpdateNextRep;

    // Waiting for the command to run the next simulation
    If (not RunRepsContinuosly) then
     begin
      RunNextRepEvent.ResetEvent;
      WaitForSingleObject(RunNextRepEvent.Handle, INFINITE);
     end;
   end;

  CalcThreadFitIndex;

  // Should be disabled in case of MCMC search!
  ExportParComboEstResults;

  UpdateSimFinished;

  SimFinished.SetEvent;
  SimFinished.ResetEvent;
 end;

Procedure TSim.UpdateNextStep;
 begin
  UpdateNextStepExternalProcedure;
 end;

Procedure TSim.UpdateNextRep;
 var
  i: Integer;
 begin
  Try
    If TTextRec(FileOutReplicate^).Mode <> 0 then
     begin
      EnterCriticalSection(FileOutReplicateWrite^);

      WriteLn(FileOutReplicate^);
      Write(FileOutReplicate^,
            LstRandSeed:0, #9,
            InitCell:0, #9);

      for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
        Write(FileOutReplicate^,
            ModelParameters.CarryingCapParams[i]:0:10, #9);

      Write(FileOutReplicate^,
            FitIndex:0:5, #9,
            FitStatsRep.nLangs:0:0, #9,
            FitStatsRep.RegResAvgRich.r2:0:5, #9,
            FitStatsRep.RegResAvgRS.r2:0:5, #9,
            FitStatsRep.KS:0:5);

      LeaveCriticalSection(FileOutReplicateWrite^);
     end;
  Except
    LeaveCriticalSection(FileOutReplicateWrite^);
   End;
 end;

Procedure TSim.UpdateSimFinished;
 begin
  // Empty function
 end;

Procedure TSim.RunNextStep;
 var
  TmpNeighborCells: TPMapCellVec;
  TmpNeighborCells2: TPMapCellVec;
  TmpNeighborLangCells: TPLangMapCellVec;

  TmpGroupCells: TPMapCellVec;
  PrevLangPatches: Array of TLangPatchInfo;
  NewLangPatches: Array of TLangPatchInfo;
  SplitLanguage: TLangPatchInfo;

  KVec: TIntVector; SumK: Int64;
  NVec: TIntVector; SumN, SumN2: Int64;
  DifNKVec: TIntVector; SumDifNK: Int64;
  IncVec: TIntVector;
  PropIncVec: TDblVector;
  IncN, DeltaPop: Int64;
  LangSeq2: TIntVector;
  CellSeq2: TIntVector;

  PopCap: Integer;
  ProbSpeciate: Double;

  TmpDbl: Double;
  TmpInt: Integer;
  TmpIntVec: TIntVector;
  TmpDblMat: TDblMatrix;

  Cell1, Cell2: PMapCell;
  D1, D2: Double;

  Found1: Boolean;
  Found2: Boolean;

  RndLangSeq: Boolean;

  PrevPatchControl: TPatchControl;
  NewPatchControl: TPatchControl;

  Language: PLanguage;
  NewLanguage: PLanguage;

  AncestorColor: TColor;

 var
  c, l, s: Integer;
  c1, c2, c3, lt, s1: Integer;

  i, j, k, w: Integer;
  p, n: Integer;
  a1, a2: Integer;

  TmpFile: TextFile;
 Label
  Lbl1;
 begin

  Try

    // Do we want to shuffle the order of languages at each time step?
    RndLangSeq:= False;

    // Re-shuffle the order of languages
    If RndLangSeq then
      LangSeq:= Sample(WrkLangs);

    // Next time step in the simulation
    TimeStep:= TimeStep + 1;

    If TimeStep > 10000 then
     begin
{$IfNDef CONSOLE}
    If AVIRec <> Nil then
      AVIRec.StopRecording;
{$EndIf}

{$IfDef CONSOLE}
      WriteLn('');
      Write('Trapped! RandomSeed: ', Self.LstRandSeed:0, '   ');
      for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
        Write(ModelParameters.CarryingCapParams[i]:0:7, '   ');
{$EndIf}

      Randomize;
      ModelParameters.RandomSeed:= RandSeed;

      ResetModel;
      Exit;

{$IfDef CONSOLE}
      WriteLn('');
      Write('Trapped! RandomSeed: ', Self.LstRandSeed:0, '   ');
      for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
        Write(ModelParameters.CarryingCapParams[i]:0:7, '   ');
      ReadLn(TimeStep);
{$EndIf}
     end;








    lt:= -1;
    While lt+1 < Length(WrkLangs) do
     begin
      lt:= lt + 1;

      Langs[WrkLangs[lt]].ChangedLastStep:= Langs[WrkLangs[lt]].ChangedThisStep;
      Langs[WrkLangs[lt]].ChangedThisStep:= False;
     end;




    // Shuffle cell order
    CellSeq:= Sample(WrkCells);
    For c2:= 0 to Length(CellSeq)-1 do
     begin
      c3:= CellSeq[c2];

      // Is the cell occupied by a language?
      If MapCells[c3].Languages = Nil then
        Continue;

      If Length(MapCells[c3].Languages) > 1 then
        c3:= c3;

      // This cell and all neighbors are already saturated
      If (MapCells[c3].PopSaturated) or
         (MapCells[c3].TotPopSize = 0) then
        Continue;

      // Copy the sequence of languages in cell c3 to vector
      SetLength(LangSeq2, Length(MapCells[c3].Languages));
      For lt:= 0 to Length(LangSeq2)-1 do
        LangSeq2[lt]:= lt;

      // Shuffle sequence of languages existing in cell c3
      If Length(LangSeq2) > 1 then
        LangSeq2:= Sample(LangSeq2);

      // For each language in cell c3, in a random order
      For lt:= 0 to Length(LangSeq2)-1 do
       begin
        // Get language id
        l:= MapCells[c3].Languages[LangSeq2[lt]].IdLang;

        // Search for the record of the cell in the list of cells occupied by the language
        c:= -1;
        For c1:= 0 to Length(Langs[l].OccupCells)-1 do
         begin
          If MapCells[c3].IdCell = Langs[l].OccupCells[c1].IdCell then
           begin
            c:= c1;
            Break;
           end;
         end;

        // Get information of sourounding cells, given a language

        // First in the vector is the current cell
        TmpNeighborCells:= GetAllNeighbors(@Langs[l].OccupCells[c]);

        // Get all neighboring cells. Cells not occupied by Langs[l] are Nil
        TmpNeighborLangCells:= GetNeighborsSameLanguage(@Langs[l], TmpNeighborCells);

        // Stores Carrying Capacity
        SumK:= 0;
        SumN:= 0;
        SetLength(KVec, Length(TmpNeighborCells));
        For i:= 0 to Length(TmpNeighborCells)-1 do
         begin
          KVec[i]:= 0;

          // Calculates carrying capacity given the environment
          // Each language may respond differently to the environment (subsistence mode)
          KVec[i]:= CalcCarryingCapacity(@Langs[l], TmpNeighborCells[i]);

          // Regional carrying capacity
          SumK:= SumK + KVec[i];

          // Total current population in the neighborhood
          SumN:= SumN + TmpNeighborCells[i].TotPopSize;
         end;

        // This cell and its neighbours are fully saturated
        If SumN = SumK then
         begin
          // Remove cell c3 from the list of cells to be checked.
          // The cell is saturated and does not need re-checking
          RemoveFromWrkCells(MapCells[c3]);

          // Set as saturated.
          MapCells[c3].PopSaturated:= True;

          // Skip loop
          Continue;;
         end;

        SumN:= 0;
        SumN2:= 0;
        SumDifNK:= 0;

        SetLength(NVec, Length(TmpNeighborLangCells));
        SetLength(DifNKVec, Length(TmpNeighborLangCells));
        SetLength(CellSeq2, Length(TmpNeighborLangCells));

        // Population size change
        SetLength(IncVec, Length(TmpNeighborCells));
        For i:= 0 to Length(TmpNeighborLangCells)-1 do
         begin
          IncVec[i]:= 0;
          NVec[i]:= 0;

          CellSeq2[i]:= i;

          If TmpNeighborLangCells[i] <> Nil then
           begin
            // Get current population size
            NVec[i]:= TmpNeighborLangCells[i].LangPopSizeInCell;

            // if it is the focal cell
            If i = 0 then
             begin
              // If the cell population is above carrying capacity
              If NVec[i] > KVec[i] then
                // Reduce population to carrying capacity
                IncVec[i]:= KVec[i] - NVec[i];
             end;

            // Regional population size
            SumN:= SumN + NVec[i] + IncVec[i];

            DeltaPop:= TmpNeighborLangCells[i].DeltaLangPop;
           end
          Else
           begin
            NVec[i]:= 0;
            DeltaPop:= 0;
           end;

          // Difference between carrying capacity and current population,
          // also take into account population change in the next cycle
          DifNKVec[i]:= KVec[i] - (NVec[i] + DeltaPop);
          SumDifNK:= SumDifNK + DifNKVec[i];

          SumN2:= SumN2 + NVec[i] + DeltaPop;
         end;

        // Population increase
        If SumDifNK > 0 then
         begin
          If ((SumN2 - SumK) <> 0) and (SumK <> 0) then
           begin
            // Regional population increase. Current population limits pop increase
            //IncN:= Trunc((r*SumN)*(1-(SumN/SumK)));

            // Regional population increase. Current and future population limits pop increase
            //IncN:= Trunc((r*SumN) * (1-(SumN2/SumK)));

            // Regional population increase. Current and future population limits pop increase
            IncN:= Round((ModelParameters.r*Langs[l].OccupCells[c].LangPopSizeInCell) * (1-(SumN2/SumK)));

            // Fix possible rounding errors at very low values
            If IncN = 0 then
              IncN:= 1;
           end
          Else
            IncN:= 0;



          // Fix any rounding error
          If IncN > SumDifNK then
            IncN:= SumDifNK;

          SumN:= 0;
          SumN2:= 0;

          SetLength(PropIncVec, Length(DifNKVec));
          While SumN < IncN do
           begin
            CellSeq2:= Sample(CellSeq2);
            For i:= 0 to Length(CellSeq2)-1 do
             begin
              PropIncVec[CellSeq2[i]]:= 0;

              // Define population change in next cycle
              If IncVec[CellSeq2[i]] < 0 then
               begin
                CellChangeDeltaPop(@Langs[l], TmpNeighborCells[CellSeq2[i]], IncVec[CellSeq2[i]]);
                IncVec[CellSeq2[i]]:= 0;
               end Else

              If DifNKVec[CellSeq2[i]] = 0 then
               begin
                Continue;
               end

              Else
               begin
                // Proportion of increase in each cell
                PropIncVec[CellSeq2[i]]:= (DifNKVec[CellSeq2[i]] / SumDifNK);

                // Actual increase in each cell
                IncVec[CellSeq2[i]]:= IncVec[CellSeq2[i]] + Round(PropIncVec[CellSeq2[i]] * IncN);

                SumN:= SumN + IncVec[CellSeq2[i]];

                If SumN > IncN then
                 begin
                  IncVec[CellSeq2[i]]:= IncVec[CellSeq2[i]] - (SumN - IncN);
                  SumN:= IncN;
                 end;

                // Define population change in next cycle
                If IncVec[CellSeq2[i]] <> 0 then
                  CellChangeDeltaPop(@Langs[l], TmpNeighborCells[CellSeq2[i]], IncVec[CellSeq2[i]]);
               end;
             end;

            // There has been no increase in population since last loop
            // Possible cause is rounding errors
            If SumN = SumN2 then
              Break;

            SumN2:= SumN;
           end;

          // Due to rounding errors, population has not increased properly
          If SumN <> IncN then
           begin
            // Distribute the new population randomly among cells
            CellSeq2:= Sample(CellSeq2);
            For i:= 0 to Length(CellSeq2)-1 do
             begin
              If PropIncVec[CellSeq2[i]] = 0 then
                Continue;
              IncVec[CellSeq2[i]]:= IncVec[CellSeq2[i]] + 1;
              SumN:= SumN + 1;
              CellChangeDeltaPop(@Langs[l], TmpNeighborCells[CellSeq2[i]], IncVec[CellSeq2[i]]);
              If SumN = IncN then
                Break;
             end;
           end;

         end Else

        // Population decrease
        If SumDifNK < 0 then
         begin
          For i:= 0 to Length(IncVec)-1 do
           begin
            If IncVec[i] <> 0 then
              CellChangeDeltaPop(@Langs[l], TmpNeighborCells[i], IncVec[i]);
           end;
         end;
       end;
     end;





    // Check if the simulation has reached a halt
    // The simulation reaches an end when all cells have been colonized by at least one language
    // So that no language can expand in any direction
    l:= 0;
    lt:= -1;
    // For each working language
    While lt+1 < Length(WrkLangs) do
     begin
      lt:= lt + 1;

      // If at least one language has been changed during this time step
      If Langs[WrkLangs[lt]].ChangedThisStep then
       begin
        // Sum one here
        l:= 1;

        // Skip the search
        Break;
       end;
     end;






    // No languages have been changed during this time step
    // This means there is no longer any dynamic in the simulation
    If l = 0 then
     begin

      // Check all existing languages to see if ANY of them have an empty cell around its range
      // If there is any cell with a sourouding empty cell, a new language emerges at this cell

      LangSeq2:= Nil;
      For i:= 0 to Length(Langs)-1 do
       begin
        If Langs[i].TotPopSize > 0 then
         begin
          SetLength(LangSeq2, Length(LangSeq2)+1);
          LangSeq2[Length(LangSeq2)-1]:= i;
         end;
       end;

      LangSeq:= Sample(LangSeq2);

      For lt:= 0 to Length(LangSeq)-1 do
       begin
        Language:= @Langs[LangSeq[lt]];

        TmpNeighborCells2:= Nil;

        // For each cell occupied by the language
        For c:= 0 to Length(Language.OccupCells)-1 do
         begin

          // Is it a cell at the border?
          If Language.OccupCells[c].BorderCell then
           begin
            // Get all neighboring cells of the border cell
            TmpNeighborCells:= Nil;
            TmpNeighborCells:= GetAllNeighbors(@Language.OccupCells[c]);

            If TmpNeighborCells = Nil then
              Continue;

            For i:= 0 to Length(TmpNeighborCells)-1 do
             begin
              // Is this neighboring cell unoccupied?
              If TmpNeighborCells[i].Languages = Nil then
               begin
                // For each already recorded cell
                For j:= 0 to Length(TmpNeighborCells2)-1 do
                 begin
                  // Has this cell already recorded to be sampled?
                  If TmpNeighborCells[i] = TmpNeighborCells2[j] then
                   begin
                    // If yes, then erase it
                    TmpNeighborCells[i]:= Nil;
                    Break;
                   end;
                 end;

                // If this possible cell has not been erased, record it
                If TmpNeighborCells[i] <> Nil then
                 begin
                  SetLength(TmpNeighborCells2, Length(TmpNeighborCells2)+1);
                  TmpNeighborCells2[Length(TmpNeighborCells2)-1]:= TmpNeighborCells[i];
                 end;
               end;
             end;
           end;
         end;

        // Is there an empty neighbor?

        // If this is true multiple new languages will be generated at each time step, as long as there are
        // enough empty neighbor cells. If this is false, then only a single language will be generated at the
        // edge of the geographic distribution
        If TmpNeighborCells2 <> Nil then
         begin;

          // Sample one of the cells occupied by the language with an empty neighbor cell
          Cell1:= TmpNeighborCells2[IntRND(0,Length(TmpNeighborCells2)-1)];

          // Create a new language with the population and cells of patch n
          NewLanguage:= CreateNewLanguage(Language,
                                          Cell1,
                                          Language.SubsistModeProp);

          NewLanguage.RecentLang:= True;

          NewLanguage.ChangedThisStep:= True;

          For i:= 0 to Length(Language.OccupCells)-1 do
           begin
            Language.OccupCells[i].MapCell.PopSaturated:= True;
            RemoveFromWrkCells(Language.OccupCells[i].MapCell^);
           end;

          SetLength(WrkCells, Length(WrkCells)+1);
          WrkCells[Length(WrkCells)-1]:= Cell1.IdCell;

          Language.ChangedThisStep:= True;

          Break;
         end

        // Impossible to create a new language, as the old language is sourounded
        Else
         begin
          // Has all languages been tested, and none of them have surounding empty cells?
          If lt = Length(LangSeq)-1 then
            // Finish the simulation!
            Halted:= True;
         end;
       end;
     end;




    // Check all languages again
    // For each existing language
    lt:= -1;
    While lt+1 < Length(WrkLangs) do
     begin
      lt:= lt + 1;

      // Progress with the next language sequentially
      l:= lt;

      // Pick a language in a random order
      l:= WrkLangs[lt];
      If RndLangSeq then
        l:= LangSeq[lt];

      // Language is good?
      If (Langs[l].NeverExisted) or (Langs[l].ExtinctLang) or
         (Langs[l].RecentLang) or (not Langs[l].ExtantLang) or
         (Langs[l].TotPopSize = 0) then
        Continue;

      // Any change in this language?
      If not Langs[l].ChangedThisStep then
        Continue;

      Langs[l].Area:= 0;
      Langs[l].Perimeter:= 0;

      // For each cell occupied by a language
      For c:= 0 to Length(Langs[l].OccupCells)-1 do
       begin

        If Langs[l].OccupCells[c].MapCell.PopSaturated then
          Continue;

        // Just colonized a new cell
        If (Langs[l].OccupCells[c].LangPopSizeInCell = 0) and
           (Langs[l].OccupCells[c].DeltaLangPop > 0) then
         begin
          // Language colonizes a new cell.

          Found1:= False;
          For i:= 0 to Length(Langs[l].OccupCells[c].MapCell.Languages)-1 do
           begin
            If Langs[l].OccupCells[c].MapCell.Languages[i].IdLang = Langs[l].IdLang then
             begin
              // Language already registered as present in that cell
              Found1:= True;
              Break;
             end;
           end;

          If Not Found1 then
           begin
            // Language not registered as present in the cell
            SetLength(Langs[l].OccupCells[c].MapCell.Languages, Length(Langs[l].OccupCells[c].MapCell.Languages)+1);
            Langs[l].OccupCells[c].MapCell.Languages[Length(Langs[l].OccupCells[c].MapCell.Languages)-1]:= @Langs[l];
           end;
         end;

        // Update language population, in each cell
        Langs[l].OccupCells[c].LangPopSizeInCell:= Langs[l].OccupCells[c].LangPopSizeInCell +
                                                    Langs[l].OccupCells[c].DeltaLangPop;
        If Langs[l].OccupCells[c].LangPopSizeInCell <= 0 then
          Langs[l].OccupCells[c].LangPopSizeInCell:= 0;

        // Update population in the cell, across languages
        Langs[l].OccupCells[c].MapCell.TotPopSize:= Langs[l].OccupCells[c].MapCell.TotPopSize +
                                                     Langs[l].OccupCells[c].DeltaLangPop;
        If Langs[l].OccupCells[c].MapCell.TotPopSize <= 0 then
          Langs[l].OccupCells[c].MapCell.TotPopSize:= 0;

        // Update language population across cells
        Langs[l].TotPopSize:= Langs[l].TotPopSize +
                              Langs[l].OccupCells[c].DeltaLangPop;
        If Langs[l].TotPopSize < 0 then
          Langs[l].TotPopSize:= 0;


        // Language is extinct at the cell
        If Langs[l].OccupCells[c].LangPopSizeInCell = 0 then
          Langs[l].OccupCells[c].IdCell:= -1;

        Langs[l].OccupCells[c].MapCell.DeltaPop:= Langs[l].OccupCells[c].MapCell.DeltaPop -
                                                   Langs[l].OccupCells[c].DeltaLangPop;

        Langs[l].OccupCells[c].DeltaLangPop:= 0;
       end;

      Langs[l].DeltaPop:= 0;

      Langs[l].Area:= Length(Langs[l].OccupCells);
      Langs[l].Perimeter:= 0;
      For c:= 0 to Langs[l].Area-1 do
       begin
        // At the edge of the continent?
        If Length(Langs[l].OccupCells[c].MapCell.Neighbors) < 6 then
          Langs[l].OccupCells[c].BorderCell:= True

        // Inland at the edge of the distribution?
        Else
          Langs[l].OccupCells[c].BorderCell:= (CountNeighborsSameLanguage(@Langs[l].OccupCells[c]) < 6);

        If Langs[l].OccupCells[c].BorderCell then
          Langs[l].Perimeter:= Langs[l].Perimeter + 1;
       end;
     end;



// Process 5: Languages have maximum population size
    lt:= -1;
    While lt+1 < Length(WrkLangs) do
     begin
      lt:= lt + 1;

      // Progress with the next language sequentially
      l:= lt;

      // Pick a language in a random order
      l:= WrkLangs[lt];
      If RndLangSeq then
        l:= LangSeq[lt];

      Language:= @Langs[l];

      // Language is good?
      If (Language.NeverExisted) or (Language.ExtinctLang) or
         (Language.RecentLang) or (not Language.ExtantLang) or
         (Language.TotPopSize = 0) then
        Continue;

      // Any change in this language?
      If not Language.ChangedThisStep then
        Continue;





      // Has the population reached its maximum population size?
      If (Language.TotPopSize > Language.MaxPopSize) and
         (Length(Language.OccupCells) > 1) then
       begin

        // If new language arise at the border
        // In the following code, below, the language stops expanding and a new language arise at the border
        TmpNeighborCells2:= Nil;

        // For each cell occupied by the language
        For c:= 0 to Length(Language.OccupCells)-1 do
         begin

          // Is it a cell at the border?
          If Language.OccupCells[c].BorderCell then
           begin
            // Get all neighboring cells of the border cell
            TmpNeighborCells:= Nil;
            TmpNeighborCells:= GetAllNeighbors(@Language.OccupCells[c]);

            If TmpNeighborCells = Nil then
              Continue;

            For i:= 0 to Length(TmpNeighborCells)-1 do
             begin
              // Is this neighboring cell unoccupied?
              If TmpNeighborCells[i].Languages = Nil then
               begin
                // For each already recorded cell
                For j:= 0 to Length(TmpNeighborCells2)-1 do
                 begin
                  // Has this cell already recorded to be sampled?
                  If TmpNeighborCells[i] = TmpNeighborCells2[j] then
                   begin
                    // If yes, then erase it
                    TmpNeighborCells[i]:= Nil;
                    Break;
                   end;
                 end;

                // If this possible cell has not been erased, record it
                If TmpNeighborCells[i] <> Nil then
                 begin
                  SetLength(TmpNeighborCells2, Length(TmpNeighborCells2)+1);
                  TmpNeighborCells2[Length(TmpNeighborCells2)-1]:= TmpNeighborCells[i];
                 end;
               end;
             end;
           end;
         end;

        // Is there an empty neighbor

        // If this is true multiple new languages will be generated at each time step, as long as there are
        // enough empty neighbor cells. If this is false, then only a single language will be generated at the
        // edge of the geographic distribution
        If TmpNeighborCells2 <> Nil then
         begin;

          // Sample one of the cells occupied by the language with an empty neighbor cell
          Cell1:= TmpNeighborCells2[IntRND(0,Length(TmpNeighborCells2)-1)];

          // Create a new language with the population and cells of patch n
          NewLanguage:= CreateNewLanguage(Language,
                                          Cell1,
                                          Language.SubsistModeProp);

          NewLanguage.RecentLang:= True;

          NewLanguage.ChangedThisStep:= True;

          For i:= 0 to Length(Language.OccupCells)-1 do
           begin
            Language.OccupCells[i].MapCell.PopSaturated:= True;
            RemoveFromWrkCells(Language.OccupCells[i].MapCell^);
           end;

          SetLength(WrkCells, Length(WrkCells)+1);
          WrkCells[Length(WrkCells)-1]:= Cell1.IdCell;

          Language.ChangedThisStep:= True;
         end

        // Impossible to create a new language, as the old language is sourounded
        Else
         begin
          // Make this false because there is no room of a new language at the edge
          Language.ChangedThisStep:= False;

          // Update the realized distribution of language population sizes
          UpdatePopSizeFrequency(Language);
         end;
//}
       end;
     end;








    For l:= 0 to Length(WrkLangs)-1 do
     begin
      If Langs[WrkLangs[l]].RecentLang then
       begin
        TotNumPatches:= TotNumPatches + Langs[WrkLangs[l]].NumPatches;
        Langs[WrkLangs[l]].RecentLang:= False;
        Langs[WrkLangs[l]].ExtantLang:= True;
       end;
     end;




    lt:= -1;
    While lt+1 < Length(WrkLangs) do
     begin
      lt:= lt + 1;

      If (not Langs[WrkLangs[lt]].ChangedLastStep) and
         (not Langs[WrkLangs[lt]].ChangedThisStep) then
       begin
        RemoveFromWrkLangs(@Langs[WrkLangs[lt]]);
        lt:= lt - 1;
        Continue;
       end;
     end;



    // All cells have been colonized at least once
    // Time to finish the simulation
    If ((TimeStep > 5000) and (Length(WrkLangs) = 0)) or (TotNumCellsColonizedOnce >= (TotNumCells-30)) then
     begin
      Halted:= True;
     end;

  Except
    On E: Exception do
     begin
      ErrMsg:= E.Message;
      {$IfNDef CONSOLE}
      //Synchronize(ShowError);
      {$EndIf}

      {$IfDef CONSOLE}
      ShowError;
      {$EndIf}
     end;
   end;
 end;

Procedure TSim.ShowError;
 begin
  If ErrMsg = '' then
   begin
    ErrorMessage('Unknown Error! Please try again.');
   end
  Else
   begin
    ErrorMessage(ErrMsg);
   end;
 end;

{$IfDef CONSOLE}
Procedure ErrorMessage(Msg: String);
 var
  i: Integer;
 begin
  WriteLn(Msg);
  ReadLn(i);
 end;
{$EndIf}







Function TSim.WritePhy(Parent: Integer; t: Integer): String;
 var
  l, z: Integer;
  fst: Boolean;
 begin
  fst:= True;
  Result:= '()';

  For l:= 0 to Length(Langs)-1 do
   begin
    If Langs[l].NeverExisted then
      Continue;

    If Langs[l].DaughterOf = Parent then
     begin

      if not fst then
        Insert(',', Result, Length(Result))
      Else
        fst:= not fst;

      If Langs[l].ExtantLang then
        z:= TimeStep - Langs[l].BornAt
      Else
        z:= Langs[l].DiedAt - Langs[l].BornAt;

      // Does species SpeciesIndex have children?
      If Langs[l].AncestorOf = Nil then
        // Because species Sp does not have any children, writes its name and close the node
        Insert('Lg' + IntToStr(Langs[l].IdLang) + ':' + IntToStr(z), Result, Length(Result))

      Else                                      // Species SpeciesIndex has children
        // Recursively call WritePhy to enter all the children of species SpeciesIndex
        Insert(WritePhy(Langs[l].IdLang, Langs[l].DiedAt) + 'Lg' + IntToStr(Langs[l].IdLang) + ':' + IntToStr(z), Result, Length(Result));
     end;
   end;
 end;


Function TSim.GetPopSize(Lang: PLanguage; MapCell: PMapCell): Integer;
 var
  i: Integer;
 begin
  Result:= 0;
  For i:= 0 to Length(MapCell.Languages)-1 do
   begin
    If MapCell.Languages[i].IdLang = Lang.IdLang then
     begin
      Result:= MapCell.Languages[i].TotPopSize;
      Break;
     end;
   end;
 end;


Function TSim.GetPopSize(LangMapCell: PLangMapCell): Integer;
 begin
  Result:= LangMapCell.LangPopSizeInCell;
 end;


Function TSim.GetDeltaPopSize(Lang: PLanguage; MapCell: PMapCell): Integer;
 var
  i: Integer;
 begin
  Result:= 0;
  For i:= 0 to Length(MapCell.Languages)-1 do
   begin
    If MapCell.Languages[i].IdLang = Lang.IdLang then
     begin
      Result:= MapCell.Languages[i].DeltaPop;
      Break;
     end;
   end;
 end;


Function TSim.GetDeltaPopSize(LangMapCell: PLangMapCell): Integer;
 begin
  Result:= LangMapCell.DeltaLangPop;
 end;

Function TSim.GetAllNeighbors(LangCell: PLangMapCell): TPMapCellVec;
 var
  i: Integer;
  TmpIntVec: TIntVector;
  TmpMapCellVec: TPMapCellVec;
 begin
  Result:= Nil;
  SetLength(Result, Length(LangCell.MapCell.Neighbors)+1);
  Result[0]:= LangCell.MapCell;

  For i:= 1 to Length(Result)-1 do
    Result[i]:= LangCell.MapCell.Neighbors[i-1];

  Exit;




  TmpMapCellVec:= Nil;
  SetLength(TmpMapCellVec, Length(Result));

  TmpIntVec:= Nil;
  SetLength(TmpIntVec, Length(Result));

  For i:= 0 to Length(Result)-1 do
   begin
    TmpMapCellVec[i]:= Result[i];
    TmpIntVec[i]:= i;
   end;

  TmpIntVec:= Sample(TmpIntVec);

  For i:= 0 to Length(Result)-1 do
   begin
    If TmpIntVec[i] > Length(Result)-1 then
      TmpIntVec[i]:= Length(Result)-1;

    Result[i]:= TmpMapCellVec[TmpIntVec[i]];

    If Result[i].IdCell = -1 then
      Result[i].IdCell:= -2;
   end;
 end;


// Given a language, occupying a certain cell, list all neighbouring cell occupied by a different language.
// Returns a vector of cells occupied by other languages
// Vector is empty (Nil) if surounding cells are occupied by the same language
Function TSim.GetNeighborsDifferentLanguage(LangCell: PLangMapCell): TPLangMapCellVec;
 var
  i, j, k: Integer;
  Lang: PLanguage;
  TmpIntVec: TIntVector;
  TmpMapCellVec: TPLangMapCellVec;
 begin
  Result:= Nil;

  Lang:= LangCell.Language;
  For i:= 0 to Length(LangCell.MapCell.Neighbors)-1 do
   begin
    For j:= 0 to Length(LangCell.MapCell.Neighbors[i].Languages)-1 do
     begin
      If LangCell.MapCell.Neighbors[i].Languages[j] <> Lang then
       begin
        For k:= 0 to Length(LangCell.MapCell.Neighbors[i].Languages[j].OccupCells)-1 do
         begin
          If LangCell.MapCell.Neighbors[i].Languages[j].OccupCells[k].IdCell =
             LangCell.MapCell.Neighbors[i].IdCell  then
           begin
            SetLength(Result, Length(Result)+1);
            Result[Length(Result)-1]:= @LangCell.MapCell.Neighbors[i].Languages[j].OccupCells[k];

            Break;
           end;
         end;

        Break;
       end;
     end;
   end;
 end;

Function TSim.CountNeighborsSameLanguage(LangCell: PLangMapCell): Integer;
 var
  i, j, k: Integer;
  Lang: PLanguage;
 begin
  Result:= 0;

  Lang:= LangCell.Language;

  For i:= 0 to Length(LangCell.MapCell.Neighbors)-1 do
   begin
    For j:= 0 to Length(LangCell.MapCell.Neighbors[i].Languages)-1 do
     begin
      If LangCell.MapCell.Neighbors[i].Languages[j] = Lang then
       begin
        For k:= 0 to Length(LangCell.MapCell.Neighbors[i].Languages[j].OccupCells)-1 do
         begin
          If LangCell.MapCell.Neighbors[i].Languages[j].OccupCells[k].IdCell =
             LangCell.MapCell.Neighbors[i].IdCell  then
           begin
            Result:= Result + 1;
            Break;
           end;
         end;

        Break;
       end;
     end;
   end;
 end;


// Get all neighboring cells. Cells not occupied by Langs[l] are Nil
Function TSim.GetNeighborsSameLanguage(Language: PLanguage; Neighbors: TPMapCellVec): TPLangMapCellVec;
 var
  i, j, k: Integer;
  TmpIntVec: TIntVector;
  TmpMapCellVec: TPLangMapCellVec;
 begin
  SetLength(Result, Length(Neighbors));

  For i:= 0 to Length(Neighbors)-1 do
   begin
    Result[i]:= Nil;
    For j:= 0 to Length(Neighbors[i].Languages)-1 do
     begin
      If Neighbors[i].Languages[j] = Language then
       begin
        For k:= 0 to Length(Neighbors[i].Languages[j].OccupCells)-1 do
         begin
          If Neighbors[i].Languages[j].OccupCells[k].IdCell =
             Neighbors[i].IdCell  then
           begin
            Result[i]:= @Neighbors[i].Languages[j].OccupCells[k];

            Break;
           end;
         end;

        Break;
       end;
     end;
   end;
 end;

// Given a language, occupying a certain cell, list all neighbouring cell occupied by the same language.
// Returns a vector of cells occupied by the same language languages
// Vector is empty (Nil) if surounding cells are occupied by different languages
Function TSim.GetNeighborsSameLanguage(LangCell: PLangMapCell): TPLangMapCellVec;
 var
  i, j, k: Integer;
  Lang: PLanguage;
  TmpIntVec: TIntVector;
  TmpMapCellVec: TPLangMapCellVec;
 begin
  SetLength(Result, 1);
  Result[0]:= LangCell;

  Lang:= LangCell.Language;

  For i:= 0 to Length(LangCell.MapCell.Neighbors)-1 do
   begin
    For j:= 0 to Length(LangCell.MapCell.Neighbors[i].Languages)-1 do
     begin
      If LangCell.MapCell.Neighbors[i].Languages[j] = Lang then
       begin
        For k:= 0 to Length(LangCell.MapCell.Neighbors[i].Languages[j].OccupCells)-1 do
         begin
          If LangCell.MapCell.Neighbors[i].Languages[j].OccupCells[k].IdCell =
             LangCell.MapCell.Neighbors[i].IdCell  then
           begin
            SetLength(Result, Length(Result)+1);
            Result[Length(Result)-1]:= @LangCell.MapCell.Neighbors[i].Languages[j].OccupCells[k];

            Break;
           end;
         end;

        Break;
       end;
     end;
   end;
 end;

Function TSim.CalcCarryingCapacity(Lang: PLanguage; MapCell: PMapCell): Integer;
 var
  i, j, k: Integer;
  Min, Max, Mid, Range, X, w: Double;
  AreaInKm2: Double;
  t: Double;
  LangMapCell: PLangMapCell;
 begin

  If MapCell.Environment <> Nil then
   begin

    // What is the environment of the cell?
    Result:= 0;

    For i:= 0 to TotNumEnvVars-1 do
      If IsNaN(MapCell.Environment[i]) then
        Exit;



    // Function between environment and carrying capacity

    If MapCell.Area = 0 then
      AreaInKm2:= 1
    Else
      AreaInKm2:= MapCell.Area;

    // Power function
    if ModelParameters.ModelFunction = TMFPower then
     begin
      // Calculate the log of expected density given the precipitation
      X:= ModelParameters.CarryingCapParams[0] + (ModelParameters.CarryingCapParams[1] * Log10(MapCell.Environment[0]));

      // Calculate the expected density in units of humans / km^2
      X:= Power(10, X);
     end Else


     //Uses and estimated density from  41. Kavanagh PH, et al. 2018 Hindcasting global population densities reveals forces enabling the origin of agriculture. Nat Hum Behav 2(7):478–484.
    if ModelParameters.ModelFunction = EstimatedDensity then
     begin
       X:= MapCell.Environment[0];
     end Else


    // Logistic function
    if ModelParameters.ModelFunction = TMFLogistic then
     begin
      X:= (ModelParameters.CarryingCapParams[2]) / (1 + EXP((-1 * Power(10, ModelParameters.CarryingCapParams[1])) * (MapCell.Environment[0] - ModelParameters.CarryingCapParams[0])));
     end Else


    // Exponential function
    if ModelParameters.ModelFunction = TMFExponential then
     begin
      X:= ModelParameters.CarryingCapParams[0] + Exp(Power(10, ModelParameters.CarryingCapParams[1]) * MapCell.Environment[0]);
      if X < 0 then
        X:= 0;
     end;


     // Kavanagh PH, et al. 2018 Hindcasting global population densities reveals forces enabling the origin of agriculture. Nat Hum Behav 2(7):478–484.
    X:= MapCell.Environment[0];

    // Calculate the expected number of humans given the size of the cell
    Result:= Round(X * AreaInKm2);

    //Almost never happens but if the result is < than 2 people, then set the
    //carrying capacity to 2
    If Result < 2 then
      Result:= 2;


    //!!!!! NOT USED IN THIS VERSION
    // Other people in the cell?
    For i:= 0 to Length(MapCell.Languages)-1 do
     begin
      If MapCell.Languages[i].IdLang = Lang.IdLang then
        Continue;

      // Maximum Interference competition !!!
      // There can be no two languages in a single cell
      Result:= 0;
      Break;
      end;

   end
  Else
    Result:= 0;
 end;

Function TSim.CalcCarryingCapacity(LangMapCell: PLangMapCell): Integer;
 begin
  Result:= CalcCarryingCapacity(LangMapCell.Language, LangMapCell.MapCell);
 end;

Function TSim.CalcCarryingCapacity(MapCell: PMapCell): Integer;
 begin
  Result:= CalcCarryingCapacity(Nil, MapCell);
 end;

Procedure TSim.UpdatePopSizeFrequency(Language: PLanguage);
 var
  BinObs, BinSim: Integer;
  i, t: Integer;
 begin
  Try
    If Language.PopSizeAccounted then
      Exit;



    // For each population size bin
    For i:= 0 to Length(MaxPopDistrib)-1 do
     begin
      // Has the language SAMPLED MAXIMUM population size fallen within this bin?
      If (Language.MaxPopSize > MaxPopDistrib[i].Min) and (Language.MaxPopSize <= MaxPopDistrib[i].Max) then
       begin
        // If so, then one language is added to the bin, therefore increasing its frequency
        MaxPopDistrib[i].SimCount:= MaxPopDistrib[i].SimCount + 1;
        Break;
       end;
     end;



    // For each population size bin
    For i:= 0 to Length(MaxPopDistrib)-1 do
     begin
      // Has the language SIMULATED population size fallen within this bin?
      If (Language.TotPopSize > MaxPopDistrib[i].Min) and (Language.TotPopSize <= MaxPopDistrib[i].Max) then
       begin
        // If so, then one language is subtracted from the bin, therefore reducing its frequency
        MaxPopDistrib[i].SimCount:= MaxPopDistrib[i].SimCount - 1;
        Break;
       end;
     end;




    // Now, recalculates the frequency distributions
    t:= 0;
    // For each bin
    For i:= 0 to Length(MaxPopDistrib)-1 do
      // Counts the number of languages inside.
      // The sum over all bins should be constant along the simulation
      t:= t + MaxPopDistrib[i].SimCount;



    // Re-calculates the frequency, which will be used as a probability
    For i:= 0 to Length(MaxPopDistrib)-1 do
     begin
      // Ratio between count and total
      MaxPopDistrib[i].SimFrequency:= MaxPopDistrib[i].SimCount / t;

      // Because bins are being subtracted, they may become negative probabilities
      // In this case, zero them
      If MaxPopDistrib[i].SimFrequency < 0 then
        MaxPopDistrib[i].SimFrequency:= 0;
     end;





    Language.PopSizeAccounted:= True;
  Except
    {$IfDef CONSOLE}
      WriteLn('Problem UpdatePopSizeFrequency');
      WriteLn(LstRandSeed, #9, TimeStep);
      WriteLn(i, #9, t, #9, MaxPopDistrib[i].ObsCount, #9, MaxPopDistrib[i].ObsFrequency);
      WriteLn(MaxPopDistrib[i].SimCount, #9, MaxPopDistrib[i].SimFrequency);
      ReadLn(i);
    {$EndIf}
   End;
 end;

Function TSim.SampleMaxPopSize: Integer;
 var
  i, j, Bin: Integer;
  W, C: TDblVector;
  SW, t: Double;
 begin
  Try
    SetLength(W, Length(MaxPopDistrib));

    SW:= 0;
    For i:= 0 to Length(W)-1 do
     begin
      W[i]:= MaxPopDistrib[i].SimFrequency;
      SW:= SW + W[i];
     end;

    SetLength(C, Length(MaxPopDistrib)+1);
    C[0]:= 0;
    For i:= 1 to Length(W) do
      C[i]:= C[i-1] + (W[i-1] / SW);

    t:= Random;

    For j:= 0 to Length(C)-2 do
     begin
      If (t >= C[j]) and (t <= C[j+1]) then
       begin
        Bin:= j;
        Break;
       end;
     end;

    Result:= Trunc(MaxPopDistrib[Bin].MaxPopVec[IntRND(0,Length(MaxPopDistrib[Bin].MaxPopVec)-1)]);
  Except
    {$IfDef CONSOLE}
      WriteLn('Problem SampleMaxPopSize');
      WriteLn(LstRandSeed, #9, TimeStep);
      WriteLn(i, #9, j, #9, t, #9, Bin);
      ReadLn(i);
    {$EndIf}
   End;
 end;

Function TSim.GetLangMapCell(Lang: PLanguage; MapCell: PMapCell): PLangMapCell;
 var
  c, z: Integer;
 begin
  Result:= Nil;

  z:= -1;
  For c:= 0 to Length(Lang.OccupCells)-1 do
   begin
    If Lang.OccupCells[c].IdCell = MapCell.IdCell then
     begin
      z:= c;
      Break;
     end;
   end;

  If z = -1 then
   begin
    c:= Length(Lang.OccupCells);
    SetLength(Lang.OccupCells, c+1);

    z:= c;

    Lang.OccupCells[z].IdCell:= MapCell.IdCell;
    Lang.OccupCells[z].PosVec:= z;

    Lang.OccupCells[z].MapCell:= MapCell;
    Lang.OccupCells[z].Language:= Lang;

    Lang.OccupCells[z].DeltaLangPop:= 0;
    Lang.OccupCells[z].K:= 0;
    Lang.OccupCells[z].LangPopSizeInCell:= 0;

    Lang.OccupCells[z].PatchTag:= -1;
    Lang.OccupCells[z].PrevPatchTag:= -1;
   end;

  Result:= @Lang.OccupCells[z];
 end;

Function TSim.GetLangMapCellVec(Lang: PLanguage; MapCellVec: TPMapCellVec): TPLangMapCellVec;
 var
  i: Integer;
 begin
  SetLength(Result, Length(MapCellVec));
  For i:= 0 to Length(Result)-1 do
   begin
    Result[i]:= GetLangMapCell(Lang, MapCellVec[i]);
   end;
 end;


Function TSim.CellChangePop(Lang: PLanguage; MapCell: PMapCell; DeltaNIndiv: Integer): PLangMapCell;
 var
  i: Integer;
  FoundLang: Boolean;
 begin
  Result:= GetLangMapCell(Lang, MapCell);

  If (MapCell.DeltaPop = 0) and
     (MapCell.TotPopSize = 0) then
   begin
    SetLength(WrkCells, Length(WrkCells)+1);
    WrkCells[Length(WrkCells)-1]:= MapCell.IdCell;
   end;

  FoundLang:= False;
  For i:= 0 to Length(MapCell.Languages)-1 do
   begin
    If MapCell.Languages[i].IdLang = Lang.IdLang then
     begin
      FoundLang:= True;
      Break;
     end;
   end;
  If Not FoundLang then
   begin
    SetLength(MapCell.Languages, Length(MapCell.Languages)+1);
    MapCell.Languages[Length(MapCell.Languages)-1]:= Lang;
   end;

  If MapCell.ColonizedOnce = False then
   begin
    MapCell.ColonizedOnce:= True;
    TotNumCellsColonizedOnce:= TotNumCellsColonizedOnce + 1;
   end;

  Lang.TotPopSize:= Lang.TotPopSize + DeltaNIndiv;
  Result.LangPopSizeInCell:= Result.LangPopSizeInCell + DeltaNIndiv;
  Result.MapCell.TotPopSize:= Result.MapCell.TotPopSize + DeltaNIndiv;
 end;

  // This function is not being used
Function TSim.CellChangePop(LangMapCell: PLangMapCell; DeltaNIndiv: Integer): PLangMapCell;
 begin
  Result:= LangMapCell;
  Result.Language.TotPopSize:= Result.Language.TotPopSize + DeltaNIndiv;
  Result.LangPopSizeInCell:= Result.LangPopSizeInCell + DeltaNIndiv;
  Result.MapCell.TotPopSize:= Result.MapCell.TotPopSize + DeltaNIndiv;
 end;

Function TSim.CellChangeDeltaPop(Lang: PLanguage; MapCell: PMapCell; DeltaNIndiv: Integer): PLangMapCell;
 var
  i: Integer;
  FoundLang: Boolean;
 begin
  If (MapCell.DeltaPop = 0) and
     (MapCell.TotPopSize = 0) then
   begin
    SetLength(WrkCells, Length(WrkCells)+1);
    WrkCells[Length(WrkCells)-1]:= MapCell.IdCell;
   end;

  Result:= GetLangMapCell(Lang, MapCell);

  FoundLang:= False;
  For i:= 0 to Length(MapCell.Languages)-1 do
   begin
    If MapCell.Languages[i].IdLang = Lang.IdLang then
     begin
      FoundLang:= True;
      Break;
     end;
   end;

  If Not FoundLang then
   begin

    If Length(MapCell.Languages) > 0 then
     begin
      {$IfDef CONSOLE}
       WriteLn('');
       WriteLn('Error! Two languages are trying to colonize the same cell...');
       WriteLn('Error! Check SumK integer overflow!');
       WriteLn('');
       WriteLn(LstRandSeed, #9, TimeStep, #9, MapCell.IdCell, #9, Lang.IdLang);
       ReadLn(i);
      {$EndIf}
     end;

    SetLength(MapCell.Languages, Length(MapCell.Languages)+1);
    MapCell.Languages[Length(MapCell.Languages)-1]:= Lang;
   end;

  If MapCell.ColonizedOnce = False then
   begin
    MapCell.ColonizedOnce:= True;
    TotNumCellsColonizedOnce:= TotNumCellsColonizedOnce + 1;
   end;

  Lang.ChangedThisStep:= True;

  Lang.DeltaPop:= Lang.DeltaPop + DeltaNIndiv;
  If Lang.DeltaPop + Lang.TotPopSize < 0 then
    Lang.DeltaPop:= -Lang.TotPopSize;

  MapCell.DeltaPop:= MapCell.DeltaPop + DeltaNIndiv;
  If MapCell.DeltaPop + MapCell.TotPopSize < 0 then
    MapCell.DeltaPop:= -MapCell.TotPopSize;

  Result.DeltaLangPop:= Result.DeltaLangPop + DeltaNIndiv;
  If Result.DeltaLangPop + Result.LangPopSizeInCell < 0 then
    Result.DeltaLangPop:= -Result.LangPopSizeInCell;
 end;

Procedure TSim.RemoveFromWrkCells(MapCell: TMapCell);
 var
  i, j: Integer;
 begin
  For i:= 0 to Length(WrkCells)-1 do
   begin
    If WrkCells[i] = MapCell.IdCell then
     begin
      For j:= i to Length(WrkCells)-2 do
       begin
        WrkCells[j]:= WrkCells[j+1];
       end;
      SetLength(WrkCells, Length(WrkCells)-1);
      Break;
     end;
   end;
 end;

Procedure TSim.RemoveFromWrkLangs(Language: PLanguage);
 var
  i, j: Integer;
 begin
  For i:= 0 to Length(WrkLangs)-1 do
   begin
    If WrkLangs[i] = Language.IdLang then
     begin
      For j:= i to Length(WrkLangs)-2 do
       begin
        WrkLangs[j]:= WrkLangs[j+1];
       end;
      SetLength(WrkLangs, Length(WrkLangs)-1);
      Break;
     end;
   end;

  UpdatePopSizeFrequency(Language);
 end;

Procedure TSim.AddtoWrkCells(MapCell: TMapCell);
 var
  i, j: Integer;
 begin
  For i:= 0 to Length(WrkCells)-1 do
    If WrkCells[i] = MapCell.IdCell then
      Exit;

  SetLength(WrkCells, Length(WrkCells)+1);
  WrkCells[Length(WrkCells)-1]:= MapCell.IdCell;
 end;

Procedure TSim.AddToWrkLangs(Language: PLanguage);
 var
  i, j: Integer;
 begin
  For i:= 0 to Length(WrkLangs)-1 do
    If WrkLangs[i] = Language.IdLang then
      Exit;

  SetLength(WrkLangs, Length(WrkLangs)+1);
  WrkLangs[Length(WrkLangs)-1]:= Language.IdLang;
 end;

// This function is not being used
Function TSim.CellChangeDeltaPop(LangMapCell: PLangMapCell; DeltaNIndiv: Integer): PLangMapCell;
 var
  i: Integer;
  FoundLang: Boolean;
 begin
  Result:= LangMapCell;

  FoundLang:= False;
  For i:= 0 to Length(LangMapCell.MapCell.Languages)-1 do
   begin
    If LangMapCell.MapCell.Languages[i].IdLang = Result.Language.IdLang then
     begin
      FoundLang:= True;
      Break;
     end;
   end;
  If Not FoundLang then
   begin
    SetLength(LangMapCell.MapCell.Languages, Length(LangMapCell.MapCell.Languages)+1);
    LangMapCell.MapCell.Languages[Length(LangMapCell.MapCell.Languages)-1]:= Result.Language;
   end;

  If LangMapCell.MapCell.ColonizedOnce = False then
   begin
    LangMapCell.MapCell.ColonizedOnce:= True;
    TotNumCellsColonizedOnce:= TotNumCellsColonizedOnce + 1;
   end;

  Result.Language.DeltaPop:= Result.Language.DeltaPop + DeltaNIndiv;
  Result.DeltaLangPop:= Result.DeltaLangPop + DeltaNIndiv;
  Result.MapCell.DeltaPop:= Result.MapCell.DeltaPop + DeltaNIndiv;
 end;

Function TSim.ChangeLanguageIdentity(AncestorLang: PLanguage): PLanguage;
 var
  i, j: Integer;
 begin
  TotNumLangs:= TotNumLangs; // Does not change
  TotNumLangsEver:= TotNumLangsEver + 1;

  If TotNumLangsEver > Length(Langs) then
    Raise Exception.Create('Maximum number of languages reached in ChangeLanguageIdentity');

  For i:= 0 to Length(WrkLangs)-1 do
   begin
    If WrkLangs[i] = AncestorLang.IdLang then
     begin
      For j:= i to Length(WrkLangs)-2 do
        WrkLangs[j]:= WrkLangs[j+1];
      Break;
     end;
   end;

  WrkLangs[Length(WrkLangs)-1]:= TotNumLangsEver-1;

  Langs[TotNumLangsEver-1]:= AncestorLang^;

  Result:= @Langs[TotNumLangsEver-1];

  SetLength(AncestorLang.AncestorOf, Length(AncestorLang.AncestorOf)+1);
  AncestorLang.AncestorOf[Length(AncestorLang.AncestorOf)-1]:= Result;

  AncestorLang.ExtantLang:= False;
  AncestorLang.TransformedLang:= True;
  AncestorLang.NeverExisted:= False;
  AncestorLang.ExtinctLang:= False;
  AncestorLang.RecentLang:= False;
  AncestorLang.DiedAt:= TimeStep;

  Result.IdLang:= TotNumLangsEver-1;

  Result.ExtantLang:= True;
  Result.TransformedLang:= False;
  Result.NeverExisted:= False;
  Result.ExtinctLang:= False;
  Result.RecentLang:= True;
  Result.BornAt:= TimeStep;
  Result.DaughterOf:= AncestorLang.IdLang;

  SetLength(Result.OccupCells, Length(Result.OccupCells));
  For i:= 0 to Length(Result.OccupCells)-1 do
   begin
    Result.OccupCells[i].Language:= Result;

    If Result.OccupCells[i].MapCell.ColonizedOnce = False then
     begin
      Result.OccupCells[i].MapCell.ColonizedOnce:= True;
      TotNumCellsColonizedOnce:= TotNumCellsColonizedOnce + 1;
     end;

    For j:= 0 to Length(Result.OccupCells[i].MapCell.Languages)-1 do
     begin
      If Result.OccupCells[i].MapCell.Languages[j] = AncestorLang then
       begin
        Result.OccupCells[i].MapCell.Languages[j]:= Result;

        Break;
       end;
     end;
   end;

  SetLength(Result.SubsistModeProp, Length(Result.SubsistModeProp));
  SetLength(Result.Patches, Length(Result.Patches));

  Result.AncestorOf:= Nil;
  Result.Color:= RGB(IntRND(0,255), IntRND(0,255), IntRND(0,255));
 end;

Function TSim.CreateNewLanguage(AncestorLang: PLanguage; SubsistModeProp: TDblVector): PLanguage;
 begin
  TotNumLangs:= TotNumLangs + 1;
  TotNumLangsEver:= TotNumLangsEver + 1;

  If TotNumLangsEver > Length(Langs) then
    Raise Exception.Create('Maximum number of languages reached in CreateNewLanguage');

  SetLength(WrkLangs, Length(WrkLangs)+1);
  WrkLangs[Length(WrkLangs)-1]:= TotNumLangsEver-1;

  Result:= @Langs[TotNumLangsEver-1];

  Result.IdLang:= TotNumLangsEver-1;
  Result.TotPopSize:= 0;
  Result.MaxPopSize:= SampleMaxPopSize;
  Result.PopSizeAccounted:= False;
  Result.DeltaPop:= 0;
  Result.NumPatches:= 0;

  Result.ExtantLang:= True;
  Result.NeverExisted:= False;
  Result.ExtinctLang:= False;
  Result.RecentLang:= False;
  Result.TransformedLang:= False;

  If DebugMode then
   begin
    If TotNumLangsEver = 1 then
      Result.Color:= RGB(255, 0, 0) else
    If TotNumLangsEver = 2 then
      Result.Color:= RGB(0, 255, 0) else
    If TotNumLangsEver = 3 then
      Result.Color:= RGB(0, 0, 255)
    else
      Result.Color:= RGB(IntRND(0,255), IntRND(0,255), IntRND(0,255));
   end
  Else
    Result.Color:= RGB(IntRND(0,255), IntRND(0,255), IntRND(0,255));

  Result.SubsistModeProp:= SubsistModeProp;
  SetLength(Result.SubsistModeProp, TotNumCarryingCapParams);

  Result.OccupCells:= Nil;
  Result.Patches:= Nil;

  Result.NumPatches:= 1;

  If AncestorLang = Nil then
    Result.DaughterOf:= -1
  Else
    Result.DaughterOf:= AncestorLang.IdLang;

  Result.ChangedThisStep:= True;
  Result.ChangedLastStep:= True;

  Result.BornAt:= TimeStep;
  Result.DiedAt:= -1;
 end;

Function TSim.CreateNewLanguage(AncestorLang: PLanguage; Cell: PMapCell; SubsistModeProp: TDblVector): PLanguage;
 var
  PopSize: Integer;
 begin
  Result:= CreateNewLanguage(AncestorLang, SubsistModeProp);

  PopSize:= CalcCarryingCapacity(Result, Cell);

  SetLength(AncestorLang.AncestorOf, Length(AncestorLang.AncestorOf)+1);
  AncestorLang.AncestorOf[Length(AncestorLang.AncestorOf)-1]:= Result;

  Result.ChangedThisStep:= True;

  Result.TotPopSize:= PopSize;
  SetLength(Result.OccupCells, 1);
  Result.OccupCells[0].IdCell:= Cell.IdCell;
  Result.OccupCells[0].PosVec:= 0;
  Result.OccupCells[0].MapCell:= Cell;
  Result.OccupCells[0].Language:= Result;
  Result.OccupCells[0].K:= 0;
  Result.OccupCells[0].LangPopSizeInCell:= PopSize;
  Result.OccupCells[0].DeltaLangPop:= 0;
  Result.OccupCells[0].PatchTag:= 0;
  Result.OccupCells[0].PrevPatchTag:= 0;
  Result.OccupCells[0].BorderCell:= True;

  SetLength(Cell.Languages, 1);
  Cell.Languages[0]:= Result;
  Cell.TotPopSize:= PopSize;
  Cell.PopSaturated:= False;

  If Cell.ColonizedOnce = False then
   begin
    Cell.ColonizedOnce:= True;
    TotNumCellsColonizedOnce:= TotNumCellsColonizedOnce + 1;
   end;

  SetLength(Result.Patches, 1);
  SetLength(Result.Patches[0].Cells, 1);
  Result.Patches[0].Cells[0]:= Cell;
 end;


Function TSim.CreateNewLanguage(AncestorLang: PLanguage; CellOccup: TPMapCellVec; SubsistModeProp: TDblVector): PLanguage;
 var
  i, j, c: Integer;
  TmpMapCellVec: TPLangMapCellVec;
 begin
  Result:= CreateNewLanguage(AncestorLang, SubsistModeProp);

  SetLength(AncestorLang.AncestorOf, Length(AncestorLang.AncestorOf)+1);
  AncestorLang.AncestorOf[Length(AncestorLang.AncestorOf)-1]:= Result;

  TmpMapCellVec:= GetLangMapCellVec(AncestorLang, CellOccup);

  SetLength(Result.OccupCells, Length(TmpMapCellVec));
  For i:= 0 to Length(TmpMapCellVec)-1 do
   begin
    Result.TotPopSize:= Result.TotPopSize + TmpMapCellVec[i].LangPopSizeInCell;

    Result.OccupCells[i].IdCell:= TmpMapCellVec[i].IdCell;
    Result.OccupCells[i].PosVec:= i;
    Result.OccupCells[i].MapCell:= TmpMapCellVec[i].MapCell;
    Result.OccupCells[i].Language:= Result;
    Result.OccupCells[i].K:= 0;
    Result.OccupCells[i].LangPopSizeInCell:= TmpMapCellVec[i].LangPopSizeInCell;
    Result.OccupCells[i].DeltaLangPop:= 0;
    Result.OccupCells[i].PatchTag:= 0;
    Result.OccupCells[i].PrevPatchTag:= 0;

    AncestorLang.TotPopSize:= AncestorLang.TotPopSize - TmpMapCellVec[i].LangPopSizeInCell;

    // Temporary flag to be "deleted" below
    TmpMapCellVec[i].IdCell:= -5;
   end;

  For i:= 0 to Length(CellOccup)-1 do
   begin
    For j:= 0 to Length(CellOccup[i].Languages)-1 do
     begin
      If CellOccup[i].Languages[j] = AncestorLang then
       begin
        CellOccup[i].Languages[j]:= Result;
        Break;
       end;
     end;

    If CellOccup[i].ColonizedOnce = False then
     begin
      CellOccup[i].ColonizedOnce:= True;
      TotNumCellsColonizedOnce:= TotNumCellsColonizedOnce + 1;
     end;
   end;


  SetLength(Result.Patches, 1);
  Result.Patches[0].Cells:= CellOccup;
  SetLength(Result.Patches[0].Cells, Length(CellOccup));

  SetLength(Result.Patches[0].PatchSeparationTime, 1);
  Result.Patches[0].PatchSeparationTime[0]:= 0;




  c:= 0;
  While c <= Length(AncestorLang.OccupCells)-1 do
   begin
    // If there is some... erase it from the language cell list
    If AncestorLang.OccupCells[c].IdCell = -5 then
     begin
      For i:= c to Length(AncestorLang.OccupCells)-2 do
       begin
        AncestorLang.OccupCells[i]:= AncestorLang.OccupCells[i+1];
       end;

      SetLength(AncestorLang.OccupCells, Length(AncestorLang.OccupCells)-1);


      If (AncestorLang.OccupCells = Nil) and (AncestorLang.IdLang <> -2) then
       begin
        TotNumLangs:= TotNumLangs - 1;
        AncestorLang.ExtantLang:= False;
        AncestorLang.ExtinctLang:= True;
        AncestorLang.DiedAt:= TimeStep;
       end;
     end

    Else
      c:= c + 1;
   end;
 end;

Function TSim.IntRND(Min, Max: Double): Integer;
 begin
  Result := Trunc(((Max - Min + 1) * (Random) ) + Min);
 end;

Function TSim.SetUpModel(ReSetMapData: Boolean = True): Boolean;
 var
  i, j, x: Integer;
  TmpEnvMat: TDblMatrix;
  TmpAreaMat: TDblMatrix;
  TmpDblMat: TDblMatrix;
  TmpStr: String;

  c, K, n: Integer;
  NewLang: PLanguage;
  NewOccupCell: PLangMapCell;
  SubsistModeProp: TDblVector;
 begin

  Result:= False;


  FitStatsRep.nLangs:= 0;

  FitStatsRep.KS:= 0;

  FitStatsRep.RegResAvgRS.b:= Nil;
  FitStatsRep.RegResAvgRS.bStd:= Nil;
  FitStatsRep.RegResAvgRS.t:= Nil;
  FitStatsRep.RegResAvgRS.SEb:= Nil;
  FitStatsRep.RegResAvgRS.Est:= Nil;
  FitStatsRep.RegResAvgRS.Res:= Nil;
  FitStatsRep.RegResAvgRS.r2:= 0;
  FitStatsRep.RegResAvgRS.r2Adj:= 0;
  FitStatsRep.RegResAvgRS.F:= 0;
  FitStatsRep.RegResAvgRS.Probr2:= 0;
  FitStatsRep.RegResAvgRS.AIC:= 0;
  FitStatsRep.RegResAvgRS.ESS:= 0;
  FitStatsRep.RegResAvgRS.KL:= 0;
  FitStatsRep.RegResAvgRS.JS:= 0;
  FitStatsRep.RegResAvgRS.Hel:= 0;
  FitStatsRep.RegResAvgRS.Euclid:= 0;
  FitStatsRep.RegResAvgRS.L1:= 0;

  FitStatsRep.RegResAvgRich.b:= Nil;
  FitStatsRep.RegResAvgRich.bStd:= Nil;
  FitStatsRep.RegResAvgRich.t:= Nil;
  FitStatsRep.RegResAvgRich.SEb:= Nil;
  FitStatsRep.RegResAvgRich.Est:= Nil;
  FitStatsRep.RegResAvgRich.Res:= Nil;
  FitStatsRep.RegResAvgRich.r2:= 0;
  FitStatsRep.RegResAvgRich.r2Adj:= 0;
  FitStatsRep.RegResAvgRich.F:= 0;
  FitStatsRep.RegResAvgRich.Probr2:= 0;
  FitStatsRep.RegResAvgRich.AIC:= 0;
  FitStatsRep.RegResAvgRich.ESS:= 0;
  FitStatsRep.RegResAvgRich.KL:= 0;
  FitStatsRep.RegResAvgRich.JS:= 0;
  FitStatsRep.RegResAvgRich.Hel:= 0;
  FitStatsRep.RegResAvgRich.Euclid:= 0;
  FitStatsRep.RegResAvgRich.L1:= 0;


  TotNumCarryingCapParams:= Length(ModelParameters.CarryingCapParams);
  If TotNumCarryingCapParams = 0 then
   begin
    ModelParameters.ModelFunction:= TMFPower;
    TotNumCarryingCapParams:= 2;
    SetLength(ModelParameters.CarryingCapParams, TotNumCarryingCapParams);
    ModelParameters.CarryingCapParams[0]:= -8.0;
    ModelParameters.CarryingCapParams[1]:= +2.6;
   end;

  SubsistModeProp:= Nil;
  SetLength(SubsistModeProp, TotNumCarryingCapParams);
  For i:= 0 to TotNumCarryingCapParams-1 do
    SubsistModeProp[i]:= 0;
  SubsistModeProp[0]:= 1;






  MaxPopVec:= OpenSimple(ModelParameters.DataPath + ModelParameters.MaxPopSizeFile);


    // Create 40 frequency bins for the distribution of language max pop size
    SetLength(MaxPopDistrib, 22);

    // Create the limits of each frequency class
    MaxPopDistrib[0].Min:= 0;
    MaxPopDistrib[0].Max:= 150;

    For i:= 1 to 17 do
      MaxPopDistrib[i].Max:= MaxPopDistrib[i-1].Max + 150;

    For i:= 18 to Length(MaxPopDistrib)-1 do
      MaxPopDistrib[i].Max:= MaxPopDistrib[i-1].Max + 500;

  For i:= 1 to Length(MaxPopDistrib)-1 do
   begin
    MaxPopDistrib[i].MaxPopVec:= Nil;
    MaxPopDistrib[i].Min:= MaxPopDistrib[i-1].Max;
   end;

  // Allocate each language size in the respective bin
  For j:= 0 to Length(MaxPopVec)-1 do
   begin
    For i:= 0 to Length(MaxPopDistrib)-1 do
     begin
      If (MaxPopVec[j,0] > MaxPopDistrib[i].Min) and (MaxPopVec[j,0] <= MaxPopDistrib[i].Max) then
       begin
        SetLength(MaxPopDistrib[i].MaxPopVec, Length(MaxPopDistrib[i].MaxPopVec)+1);
        MaxPopDistrib[i].MaxPopVec[Length(MaxPopDistrib[i].MaxPopVec)-1]:= MaxPopVec[j,0];

        Break;
       end;
     end;
   end;


  // Calculate the frequency of each bin
  For i:= 0 to Length(MaxPopDistrib)-1 do
   begin
    MaxPopDistrib[i].ObsCount:= Length(MaxPopDistrib[i].MaxPopVec);
    MaxPopDistrib[i].ObsFrequency:= MaxPopDistrib[i].ObsCount / Length(MaxPopVec);

    MaxPopDistrib[i].SimCount:= 0;
    MaxPopDistrib[i].SimFrequency:= MaxPopDistrib[i].ObsFrequency;
   end;



  If ModelParameters.AreaParamFile <> '' then
    TmpAreaMat:= OpenSimple(ModelParameters.DataPath + ModelParameters.AreaParamFile);

  TmpEnvMat:= OpenSimple(ModelParameters.DataPath + ModelParameters.EnvDataFile);
  TotNumEnvVars:= Length(TmpEnvMat[0]);

  If ReSetMapData then
   begin
    MapShapeSmallRes:= ReadShapeFile(ModelParameters.DataPath + ModelParameters.MapShapeSmallResFile);
    MapShape:= ReadShapeFile(ModelParameters.DataPath + ModelParameters.MapShapeFile);

    TotNumCells:= Length(MapShape.GEOData.ShpPolygon.Polygon);

    MapCells:= Nil;
    SetLength(MapCells, TotNumCells);

    CellSeq:= Nil;
    SetLength(CellSeq, TotNumCells);
    For i:= 0 to Length(MapCells)-1 do
     begin
      MapCells[i].IdCell:= i;

      MapCells[i].Coord:= Centroid(MapShape.GEOData.ShpPolygon.Polygon[i].Part[0].Vertices);

      MapCells[i].Environment:= TmpEnvMat[i];

      If TmpAreaMat <> Nil then
        MapCells[i].Area:= TmpAreaMat[i,0];

      SetLength(MapCells[i].Environment, Length(MapCells[i].Environment));

      CellSeq[i]:= i;
     end;

    For i:= 0 to TotNumCells-1 do
     begin
      For j:= i+1 to TotNumCells-1 do
       begin
        If Distance(MapCells[i].Coord, MapCells[j].Coord) <= ModelParameters.NeighborDistance then
         begin
          SetLength(MapCells[i].Neighbors, Length(MapCells[i].Neighbors)+1);
          MapCells[i].Neighbors[Length(MapCells[i].Neighbors)-1]:= @MapCells[j];

          SetLength(MapCells[j].Neighbors, Length(MapCells[j].Neighbors)+1);
          MapCells[j].Neighbors[Length(MapCells[j].Neighbors)-1]:= @MapCells[i];
         end;
       end;
     end;
   end;

 end;


Function TSim.ReSetModel: Boolean;
 var
  i, j, x: Integer;
  TmpEnvMat: TDblMatrix;
  TmpAreaMat: TDblMatrix;
  TmpStr: String;

  K, n: Integer;
  NewLang: PLanguage;
  NewOccupCell: PLangMapCell;
  SubsistModeProp: TDblVector;
 begin
  // Creates a new random seed
  If ModelParameters.RandomSeed = -1 then
    Randomize
  Else
    RandSeed:= ModelParameters.RandomSeed;

  // Stores the random seed at the begining of the simulation
  LstRandSeed:= RandSeed;

  Result:= False;

  Halted:= False;
  TimeStep:= 0;
  TotNumLangs:= 0;
  TotNumLangsEver:= 0;


  ReScaledRep.AvRangeSizePerCell:= Nil;
  ReScaledRep.RichPerCell:= Nil;
  ReScaledRep.RSFD:= Nil;

  FitStatsRep.RegResAvgRS.b:= Nil;
  FitStatsRep.RegResAvgRS.bStd:= Nil;
  FitStatsRep.RegResAvgRS.t:= Nil;
  FitStatsRep.RegResAvgRS.SEb:= Nil;
  FitStatsRep.RegResAvgRS.Est:= Nil;
  FitStatsRep.RegResAvgRS.Res:= Nil;

  FitStatsRep.RegResAvgRich.b:= Nil;
  FitStatsRep.RegResAvgRich.bStd:= Nil;
  FitStatsRep.RegResAvgRich.t:= Nil;
  FitStatsRep.RegResAvgRich.SEb:= Nil;
  FitStatsRep.RegResAvgRich.Est:= Nil;
  FitStatsRep.RegResAvgRich.Res:= Nil;

  FitIndex:= 0;

  j:= ModelParameters.NumInitLang;

  For i:= 0 to Length(MaxPopDistrib)-1 do
   begin
    MaxPopDistrib[i].SimCount:= MaxPopDistrib[i].ObsCount;
    MaxPopDistrib[i].SimFrequency:= MaxPopDistrib[i].ObsFrequency;
   end;

  Langs:= Nil;
  SetLength(Langs, 1000000); // maximum number of languages in a singulation

  LangSeq:= Nil;
  SetLength(LangSeq, Length(Langs));
  For i:= 0 to Length(Langs)-1 do
   begin
    Langs[i].ExtantLang:= False;
    Langs[i].NeverExisted:= True;
    Langs[i].ExtinctLang:= False;
    Langs[i].RecentLang:= False;
    Langs[i].TransformedLang:= False;
    Langs[i].IdLang:= -99999;
    Langs[i].BornAt:= -1;
    Langs[i].DaughterOf:= -99999;
    Langs[i].DiedAt:= -1;
    LangSeq[i]:= i;
   end;

  For i:= 0 to Length(MapCells)-1 do
   begin
    MapCells[i].Languages:= Nil;
    MapCells[i].TotPopSize:= 0;
    MapCells[i].DeltaPop:= 0;
    MapCells[i].PopSaturated:= False;
    MapCells[i].ColonizedOnce:= False;
   end;

  TotNumCellsColonizedOnce:= 0;

  WrkCells:= Nil;
  WrkLangs:= Nil;

  j:= 1;
  For i:= 1 to j do
   begin
    NewLang:= CreateNewLanguage(Nil, SubsistModeProp);

    // Seed Cell defined randomly
    If ModelParameters.InitCell = -1 then
      Self.InitCell:= IntRND(0, TotNumCells-1)
    Else
      Self.InitCell:= ModelParameters.InitCell;

    While True do
     begin
      If (Self.InitCell > Length(MapCells)-1) or
         (MapCells[Self.InitCell].Environment = Nil) or
         (IsNaN(MapCells[Self.InitCell].Environment[0])) then
       begin
        Self.InitCell:= IntRND(0, TotNumCells-1)
       end
      Else
        Break;
     end;

    K:= CalcCarryingCapacity(NewLang, @MapCells[Self.InitCell]);
    n:= ModelParameters.InitPopSize;
    NewOccupCell:= CellChangePop(NewLang, @MapCells[Self.InitCell], n);
    NewOccupCell.K:= K;

    NewLang.TotPopSize:= NewOccupCell.LangPopSizeInCell;

    For x:= 2 to ModelParameters.PopPerLanguage do
     begin
      // Seed Cell defined randomly
      If ModelParameters.InitCell = -1 then
        Self.InitCell:= IntRND(0, TotNumCells-1)
      Else
        Self.InitCell:= ModelParameters.InitCell;

      While True do
       begin
        If (Self.InitCell > Length(MapCells)-1) or
           (MapCells[Self.InitCell].Environment = Nil) or
           (IsNaN(MapCells[Self.InitCell].Environment[0])) then
         begin
          Self.InitCell:= IntRND(0, TotNumCells-1)
         end
        Else
          Break;
       end;

      K:= CalcCarryingCapacity(NewLang, @MapCells[Self.InitCell]);
      n:= ModelParameters.InitPopSize;
      NewOccupCell:= CellChangePop(NewLang, @MapCells[Self.InitCell], n);
      NewOccupCell.K:= K;

      NewLang.TotPopSize:= NewLang.TotPopSize + NewOccupCell.LangPopSizeInCell;
     end;
   end;

  TotNumPatches:= TotNumLangs;

  TimeStep:= 0;

  ModelSet:= True;




{$IfNDef CONSOLE}
  FrmMap:= TFrmMap.Create(FrmControl);
  FrmMap.SetUpMap(MapShape);
  FrmMap.Show;

  FrmMap.UpdateMap;

{$EndIf}






  Result:= True;
 end;










Function TSim.ReScale: TModelReScaleResults;
 var
  ObsPresAbs: TDataSet;
  NewPresAbs: TDataSet;
  c, p, s, i, j, k, l: Integer;
  F: TextFile;
 begin
  Try
    ObsPresAbs.NumValues:= Nil;
    ObsPresAbs.NumVarNames:= Nil;

    SetLength(ObsPresAbs.NumVarNames, TotNumLangs);
    SetLength(ObsPresAbs.NumValues, TotNumCells, TotNumLangs);

    i:= 0;
    For l:= 0 to Length(Langs)-1 do
     begin
      If (Langs[l].NeverExisted) or (Langs[l].ExtinctLang) or
         (Langs[l].RecentLang) or (not Langs[l].ExtantLang) then
       begin
        Continue;
       end;

      ObsPresAbs.NumVarNames[i]:= 'Lg' + IntToStr(i);

      For k:= 0 to Length(Langs[l].OccupCells)-1 do
        If (Langs[l].OccupCells[k].IdCell >= 0) and (Langs[l].OccupCells[k].IdCell < TotNumCells) then
          ObsPresAbs.NumValues[Langs[l].OccupCells[k].IdCell, i]:= 1;

      i:= i + 1;
     end;






    SetLength(Result.AvRangeSizePerCell, Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon));
    SetLength(Result.RichPerCell, Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon));
    SetLength(Result.RSFD, Length(ObsPresAbs.NumValues[0]));
    SetLength(NewPresAbs.NumValues, Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon), Length(ObsPresAbs.NumValues[0]));
    SetLength(NewPresAbs.NumVarNames, Length(ObsPresAbs.NumValues[0]));

    For i:= 0 to Length(NewPresAbs.NumValues)-1 do
     begin
      Result.AvRangeSizePerCell[i]:= 0;
      Result.RichPerCell[i]:= 0;
      For j:= 0 to Length(NewPresAbs.NumValues[i])-1 do
        NewPresAbs.NumValues[i,j]:= 0;
     end;

    For i:= 0 to Length(NewPresAbs.NumVarNames)-1 do
     begin
      NewPresAbs.NumVarNames[i]:= ObsPresAbs.NumVarNames[i];
      Result.RSFD[i]:= 0;
     end;

    For c:= 0 to Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon)-1 do
     begin
      For p:= 0 to Length(MapShape.GEOData.ShpPolygon.Polygon)-1 do
       begin
        If PointInPolygon(MapShape.GEOData.ShpPolygon.Polygon[p].Centroid,
                          MapShapeSmallRes.GEOData.ShpPolygon.Polygon[c].Part[0].Vertices) then
         begin
          For s:= 0 to Length(NewPresAbs.NumVarNames)-1 do
           begin
            If ObsPresAbs.NumValues[p,s] = 1 then
             begin
              NewPresAbs.NumValues[c,s]:= 1;
             end;
           end;
         end;
       end;
     end;




    For s:= 0 to Length(NewPresAbs.NumVarNames)-1 do
     begin
      c:= 0;
      For i:= 0 to Length(NewPresAbs.NumValues)-1 do
        c:= c + Trunc(NewPresAbs.NumValues[i,s]);

      Result.RSFD[s]:= c;
     end;



    // Calculate language richness in each cell
    For i:= 0 to Length(NewPresAbs.NumValues)-1 do
     begin
      For s:= 0 to Length(NewPresAbs.NumVarNames)-1 do
       begin
        If NewPresAbs.NumValues[i,s] > 0 then
          Result.RichPerCell[i]:= Result.RichPerCell[i] + 1;
       end;
     end;






    // Calculate mean range size in each cell
    For i:= 0 to Length(NewPresAbs.NumValues)-1 do
     begin
      c:= 0;
      For s:= 0 to Length(NewPresAbs.NumVarNames)-1 do
       begin
        c:= c + Trunc(NewPresAbs.NumValues[i,s]);
        Result.AvRangeSizePerCell[i]:= Result.AvRangeSizePerCell[i] + (NewPresAbs.NumValues[i,s] * Result.RSFD[s]);
       end;

      If c > 0 then
       begin
        Result.AvRangeSizePerCell[i]:= Result.AvRangeSizePerCell[i] / c;
       end
      Else
       begin
        Result.AvRangeSizePerCell[i]:= 0;
       end;
     end;
  Except
    s:= 0;
   End;
 end;


Procedure TSim.ClearParComboResults;
 var
  i, n: Integer;
 begin
  n:= Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon);

  SetLength(ParComboResults.SumRichPerCell, n);
  SetLength(ParComboResults.Sum2RichPerCell, n);
  SetLength(ParComboResults.SumAvRangeSizePerCell, n);
  SetLength(ParComboResults.Sum2AvRangeSizePerCell, n);

  For i:= 0 to n-1 do
   begin
    ParComboResults.SumRichPerCell[i]:= 0;
    ParComboResults.Sum2RichPerCell[i]:= 0;
    ParComboResults.SumAvRangeSizePerCell[i]:= 0;
    ParComboResults.Sum2AvRangeSizePerCell[i]:= 0;
   end;

  ParComboResults.AcumRSFD:= Nil;
  ParComboResults.NumLangs:= Nil;

  ParComboResults.RepCount:= 0;
 end;


Procedure TSim.StoreRepResults;
 var
  i, n: Integer;
 begin
  n:= Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon);

  For i:= 0 to n-1 do
   begin
    ParComboResults.SumRichPerCell[i]:= ParComboResults.SumRichPerCell[i] + ReScaledRep.RichPerCell[i];
    ParComboResults.Sum2RichPerCell[i]:= ParComboResults.Sum2RichPerCell[i] + Sqr(ReScaledRep.RichPerCell[i]);
    ParComboResults.SumAvRangeSizePerCell[i]:= ParComboResults.SumAvRangeSizePerCell[i] + ReScaledRep.AvRangeSizePerCell[i];
    ParComboResults.Sum2AvRangeSizePerCell[i]:= ParComboResults.Sum2AvRangeSizePerCell[i] + Sqr(ReScaledRep.AvRangeSizePerCell[i]);
   end;

  n:= Length(ParComboResults.AcumRSFD);
  SetLength(ParComboResults.AcumRSFD, Length(ParComboResults.AcumRSFD) + Length(ReScaledRep.RSFD));

  For i:= n to Length(ParComboResults.AcumRSFD)-1 do
    ParComboResults.AcumRSFD[i]:= ReScaledRep.RSFD[i-n];

  SetLength(ParComboResults.NumLangs, Length(ParComboResults.NumLangs)+1);
  ParComboResults.NumLangs[Length(ParComboResults.NumLangs)-1]:= TotNumLangs;

  ParComboResults.RepCount:= ParComboResults.RepCount + 1;
 end;



Function TSim.CalcRepFitStats: TModelFitStatsResults;
 begin
  Result.nLangs:= TotNumLangs;

  If ModelParameters.ObsRichMap <> Nil then
   begin
    Try
      // Regression between observed and predicted average range size in each cell
      Result.RegResAvgRich:= Regression(ModelParameters.ObsRichMap, VectorToMatrix(ReScaledRep.RichPerCell));
    Except

     End;
   end;

  If ModelParameters.ObsAvgRSMap <> Nil then
   begin
    Try
      // Regression between observed and predicted average range size in each cell
      Result.RegResAvgRS:= Regression(ModelParameters.ObsAvgRSMap, VectorToMatrix(ReScaledRep.AvRangeSizePerCell));
    Except

     End;
   end;

  // KS test between observed and predicted range size frequency distributions
  If ReScaledRep.RSFD <> Nil then
    Result.KS:= KSTest(ReScaledRep.RSFD, ModelParameters.ObsRSFD);

 end;

Procedure TSim.CalcAvgRepFitStats;
 var
  i: Integer;
 begin
  If (FitStatsRep.nLangs = 0) or
     IsNaN(FitStatsRep.RegResAvgRich.r2) then
   begin
    ActNumReps:= ActNumReps - 1;
    Exit;
   end;

  AvgFitStats.nLangs:= AvgFitStats.nLangs + FitStatsRep.nLangs;

  AvgFitStats.KS:= AvgFitStats.KS + FitStatsRep.KS;

  AvgFitStats.RegResAvgRS.r2:= AvgFitStats.RegResAvgRS.r2 + FitStatsRep.RegResAvgRS.r2;
  AvgFitStats.RegResAvgRS.r2Adj:= AvgFitStats.RegResAvgRS.r2Adj + FitStatsRep.RegResAvgRS.r2Adj;
  AvgFitStats.RegResAvgRS.F:= AvgFitStats.RegResAvgRS.F + FitStatsRep.RegResAvgRS.F;
  AvgFitStats.RegResAvgRS.Probr2:= AvgFitStats.RegResAvgRS.Probr2 + FitStatsRep.RegResAvgRS.Probr2;
  AvgFitStats.RegResAvgRS.AIC:= AvgFitStats.RegResAvgRS.AIC + FitStatsRep.RegResAvgRS.AIC;
  AvgFitStats.RegResAvgRS.ESS:= AvgFitStats.RegResAvgRS.ESS + FitStatsRep.RegResAvgRS.ESS;
  AvgFitStats.RegResAvgRS.KL:= AvgFitStats.RegResAvgRS.KL + FitStatsRep.RegResAvgRS.KL;
  AvgFitStats.RegResAvgRS.JS:= AvgFitStats.RegResAvgRS.JS + FitStatsRep.RegResAvgRS.JS;
  AvgFitStats.RegResAvgRS.Hel:= AvgFitStats.RegResAvgRS.Hel + FitStatsRep.RegResAvgRS.Hel;
  AvgFitStats.RegResAvgRS.Euclid:= AvgFitStats.RegResAvgRS.Euclid + FitStatsRep.RegResAvgRS.Euclid;
  AvgFitStats.RegResAvgRS.L1:= AvgFitStats.RegResAvgRS.L1 + FitStatsRep.RegResAvgRS.L1;

  AvgFitStats.RegResAvgRich.r2:= AvgFitStats.RegResAvgRich.r2 + FitStatsRep.RegResAvgRich.r2;
  AvgFitStats.RegResAvgRich.r2Adj:= AvgFitStats.RegResAvgRich.r2Adj + FitStatsRep.RegResAvgRich.r2Adj;
  AvgFitStats.RegResAvgRich.F:= AvgFitStats.RegResAvgRich.F + FitStatsRep.RegResAvgRich.F;
  AvgFitStats.RegResAvgRich.Probr2:= AvgFitStats.RegResAvgRich.Probr2 + FitStatsRep.RegResAvgRich.Probr2;
  AvgFitStats.RegResAvgRich.AIC:= AvgFitStats.RegResAvgRich.AIC + FitStatsRep.RegResAvgRich.AIC;
  AvgFitStats.RegResAvgRich.ESS:= AvgFitStats.RegResAvgRich.ESS + FitStatsRep.RegResAvgRich.ESS;
  AvgFitStats.RegResAvgRich.KL:= AvgFitStats.RegResAvgRich.KL + FitStatsRep.RegResAvgRich.KL;
  AvgFitStats.RegResAvgRich.JS:= AvgFitStats.RegResAvgRich.JS + FitStatsRep.RegResAvgRich.JS;
  AvgFitStats.RegResAvgRich.Hel:= AvgFitStats.RegResAvgRich.Hel + FitStatsRep.RegResAvgRich.Hel;
  AvgFitStats.RegResAvgRich.Euclid:= AvgFitStats.RegResAvgRich.Euclid + FitStatsRep.RegResAvgRich.Euclid;
  AvgFitStats.RegResAvgRich.L1:= AvgFitStats.RegResAvgRich.L1 + FitStatsRep.RegResAvgRich.L1;
 end;

Function TSim.CalcRepFitIndex: Double;
 var
  t: Double;
 begin
  If IsNan(FitStatsRep.RegResAvgRich.r2) then
   begin
    Result:= NaN;
    Exit;
   end;

  // Fit index is an arbitrary value that combines multiple measures of fit
  // Values closer to zero indicate best fit
  FitIndex:= 0;

  // r2 betwee observed and estimated language richness
  FitIndex:= FitIndex + ((1-FitStatsRep.RegResAvgRich.r2) * 0.5);

  // Number of languages
  t:= Abs(Length(ModelParameters.ObsRSFD) - FitStatsRep.nLangs);  // Absoluute difference in number of langauges
  t:= t / Length(ModelParameters.ObsRSFD);                        // Proportion between difference and number of languages

  // Calculate fit index
  FitIndex:= FitIndex + (t * 0.5);



  FitIndex:= 1 - FitIndex;
  Result:= FitIndex;
 end;

Function TSim.CalcThreadFitIndex: Double;
 var
  i, n: Integer;
  t: Double;
  v: TDblVector;
 begin
  If ActNumReps = 0 then
   begin
    Result:= 0;
    Exit;
   end;

  AvgFitStats.nLangs:= AvgFitStats.nLangs / ActNumReps;

  AvgFitStats.KS:= AvgFitStats.KS / ActNumReps;

  // The average of r2 of AvRS
  AvgFitStats.RegResAvgRS.r2:= AvgFitStats.RegResAvgRS.r2 / ActNumReps;
  AvgFitStats.RegResAvgRS.r2Adj:= AvgFitStats.RegResAvgRS.r2Adj / ActNumReps;
  AvgFitStats.RegResAvgRS.F:= AvgFitStats.RegResAvgRS.F / ActNumReps;
  AvgFitStats.RegResAvgRS.Probr2:= AvgFitStats.RegResAvgRS.Probr2 / ActNumReps;
  AvgFitStats.RegResAvgRS.AIC:= AvgFitStats.RegResAvgRS.AIC / ActNumReps;
  AvgFitStats.RegResAvgRS.ESS:= AvgFitStats.RegResAvgRS.ESS / ActNumReps;
  AvgFitStats.RegResAvgRS.KL:= AvgFitStats.RegResAvgRS.KL / ActNumReps;
  AvgFitStats.RegResAvgRS.JS:= AvgFitStats.RegResAvgRS.JS / ActNumReps;
  AvgFitStats.RegResAvgRS.Hel:= AvgFitStats.RegResAvgRS.Hel / ActNumReps;
  AvgFitStats.RegResAvgRS.Euclid:= AvgFitStats.RegResAvgRS.Euclid / ActNumReps;
  AvgFitStats.RegResAvgRS.L1:= AvgFitStats.RegResAvgRS.L1 / ActNumReps;



  // Average of r2 of the replicates
  AvgFitStats.RegResAvgRich.r2:= AvgFitStats.RegResAvgRich.r2 / ActNumReps;
  AvgFitStats.RegResAvgRich.r2Adj:= AvgFitStats.RegResAvgRich.r2Adj / ActNumReps;
  AvgFitStats.RegResAvgRich.F:= AvgFitStats.RegResAvgRich.F / ActNumReps;
  AvgFitStats.RegResAvgRich.Probr2:= AvgFitStats.RegResAvgRich.Probr2 / ActNumReps;
  AvgFitStats.RegResAvgRich.AIC:= AvgFitStats.RegResAvgRich.AIC / ActNumReps;
  AvgFitStats.RegResAvgRich.ESS:= AvgFitStats.RegResAvgRich.ESS / ActNumReps;
  AvgFitStats.RegResAvgRich.KL:= AvgFitStats.RegResAvgRich.KL / ActNumReps;
  AvgFitStats.RegResAvgRich.JS:= AvgFitStats.RegResAvgRich.JS / ActNumReps;
  AvgFitStats.RegResAvgRich.Hel:= AvgFitStats.RegResAvgRich.Hel / ActNumReps;
  AvgFitStats.RegResAvgRich.Euclid:= AvgFitStats.RegResAvgRich.Euclid / ActNumReps;
  AvgFitStats.RegResAvgRich.L1:= AvgFitStats.RegResAvgRich.L1 / ActNumReps;
  t:= AvgFitStats.RegResAvgRich.r2;



  // OR


  // Average the language richness of the replicates, and calculate one single r2
  n:= Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon);
  SetLength(v, n);

  For i:= 0 to n-1 do
    ParComboResults.SumRichPerCell[i]:= ParComboResults.SumRichPerCell[i] / ParComboResults.RepCount;

  AvgFitStats.RegResAvgRich:= Regression(ModelParameters.ObsRichMap, VectorToMatrix(ParComboResults.SumRichPerCell));
  t:= AvgFitStats.RegResAvgRich.r2;


  For i:= 0 to n-1 do
    ParComboResults.SumAvRangeSizePerCell[i]:= ParComboResults.SumAvRangeSizePerCell[i] / ParComboResults.RepCount;

  AvgFitStats.RegResAvgRS:= Regression(ModelParameters.ObsAvgRSMap, VectorToMatrix(ParComboResults.SumAvRangeSizePerCell));





  // Fit index is an arbitrary value that combines multiple measures of fit

  // Values closer to zero indicate best fit
  FitIndex:= 0;

  // r2 between observed and estimated richness map
  FitIndex:= FitIndex + ((1-t) * 0.50);





  // Number of languages
  t:= Abs(Length(ModelParameters.ObsRSFD) - AvgFitStats.nLangs);  // Absoluute difference in number of langauges
  t:= t / Length(ModelParameters.ObsRSFD);                        // Proportion between difference and number of languages



  FitIndex:= FitIndex + (t * 0.5);


  FitIndex:= 1 - FitIndex;
  Result:= FitIndex;
 end;




Procedure TSim.ExportParComboEstResults;
 var
  i, n: Integer;
  m, sd: Double;
  FileOutSim: TextFile;
 begin
  n:= Length(MapShapeSmallRes.GEOData.ShpPolygon.Polygon);

  AssignFile(FileOutSim, ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - xEstRichMap.txt');
  ReWrite(FileOutSim);
  Write(FileOutSim,
        'AvRich' + #9 +
        'StdDevRich');

  For i:= 0 to n-1 do
   begin
    WriteLn(FileOutSim);

    m:= ParComboResults.SumRichPerCell[i];
    Sd:= ((ParComboResults.Sum2RichPerCell[i]*ParComboResults.RepCount) / ParComboResults.RepCount) - Sqr(m);
    Sd:= Sqrt(Sd);

    Write(FileOutSim,
          m:0:5, #9,
          Sd:0:5);
   end;

  CloseFile(FileOutSim);







  AssignFile(FileOutSim, ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - xEstAvRSFDMap.txt');
  ReWrite(FileOutSim);
  Write(FileOutSim,
        'AvRSFD' + #9 +
        'StdDevRSFD');

  For i:= 0 to n-1 do
   begin
    WriteLn(FileOutSim);

    m:= ParComboResults.SumAvRangeSizePerCell[i];
    Sd:= ((ParComboResults.Sum2AvRangeSizePerCell[i]*ParComboResults.RepCount) / ParComboResults.RepCount) - Sqr(m);
    Sd:= Sqrt(Sd);

    Write(FileOutSim,
          m:0:5, #9,
          Sd:0:5);
   end;

  CloseFile(FileOutSim);







  AssignFile(FileOutSim, ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - xEstRSFD.txt');
  ReWrite(FileOutSim);
  Write(FileOutSim,
        'RangeSize');

  n:= Length(ParComboResults.AcumRSFD);
  For i:= 0 to n-1 do
   begin
    WriteLn(FileOutSim);

    Write(FileOutSim,
          ParComboResults.AcumRSFD[i]:0:0);
   end;

  CloseFile(FileOutSim);





  AssignFile(FileOutSim, ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - xEstN.txt');
  ReWrite(FileOutSim);
  Write(FileOutSim,
        'NumLangs');

  n:= Length(ParComboResults.NumLangs);
  For i:= 0 to n-1 do
   begin
    WriteLn(FileOutSim);

    Write(FileOutSim,
          ParComboResults.NumLangs[i]:0:0);
   end;

  CloseFile(FileOutSim);
 end;




Procedure TSim.Randomize;
 var
   Counter: Int64;
 begin
   if QueryPerformanceCounter(Counter) then
     RandSeed := Counter
   else
     RandSeed := GetTickCount;
 end;

Function TSim.Random: Extended;
 const
   two2neg32: double = ((1.0/$10000) / $10000);
 var
   Temp: Longint;
   F: Extended;
 begin
   Temp := RandSeed * $08088405 + 1;
   RandSeed := Temp;
   F := Int64(Cardinal(Temp));
   Result := F * two2neg32;
 end;

Function TSimParallel.Execute(ModelParameters: TModelParameters;
                              nReps: Integer = 120;
                              nThreads: Integer = -1): Double;
 var
  i, j, c, d: Integer;
  s: Double;
  SysInfo: SYSTEM_INFO;
  RepsPerThreads: Integer;
  TmpEnvMat: TDblMatrix;
  TmpAreaMat: TDblMatrix;
 begin
  Self.ModelParameters:= ModelParameters;

  If nThreads <= 0 then
   begin
    GetSystemInfo(SysInfo);
    nThreads:= SysInfo.dwNumberOfProcessors;
   end;

  If nThreads > nReps then
    nThreads:= nReps;

  RepsPerThreads:= nReps div nThreads;



// Open dataset only once, as it is slow to be re-opened by each thread
  If MapCells = Nil then
   begin
    If ModelParameters.AreaParamFile <> '' then
      TmpAreaMat:= OpenSimple(ModelParameters.DataPath + ModelParameters.AreaParamFile);

    TmpEnvMat:= OpenSimple(ModelParameters.DataPath + ModelParameters.EnvDataFile);
    TotNumEnvVars:= Length(TmpEnvMat[0]);

    MapShape:= ReadShapeFile(ModelParameters.DataPath + ModelParameters.MapShapeFile);
    MapShapeSmallRes:= ReadShapeFile(ModelParameters.DataPath + ModelParameters.MapShapeSmallResFile);

    TotNumCells:= Length(MapShape.GEOData.ShpPolygon.Polygon);

    MapCells:= Nil;
    SetLength(MapCells, TotNumCells);

    CellSeq:= Nil;
    SetLength(CellSeq, TotNumCells);
    For i:= 0 to Length(MapCells)-1 do
     begin
      MapCells[i].IdCell:= i;

      MapCells[i].Coord:= Centroid(MapShape.GEOData.ShpPolygon.Polygon[i].Part[0].Vertices);

      MapCells[i].Environment:= TmpEnvMat[i];

      If TmpAreaMat <> Nil then
        MapCells[i].Area:= TmpAreaMat[i,0];

      SetLength(MapCells[i].Environment, Length(MapCells[i].Environment));

      CellSeq[i]:= i;
     end;

    For i:= 0 to TotNumCells-1 do
     begin
      For j:= i+1 to TotNumCells-1 do
       begin
        If Distance(MapCells[i].Coord, MapCells[j].Coord) <= ModelParameters.NeighborDistance then
         begin
          SetLength(MapCells[i].Neighbors, Length(MapCells[i].Neighbors)+1);
          MapCells[i].Neighbors[Length(MapCells[i].Neighbors)-1]:= @MapCells[j];

          SetLength(MapCells[j].Neighbors, Length(MapCells[j].Neighbors)+1);
          MapCells[j].Neighbors[Length(MapCells[j].Neighbors)-1]:= @MapCells[i];
         end;
       end;
     end;
   end;
//



  If (ModelParameters.OutPutFolder <> '') and (ModelParameters.SimName <> '') then
   begin
    InitializeCriticalSection(FileOutReplicateWrite);

    AssignFile(FileOutReplicate, ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - Reps.txt');
    If FileExists(ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - Reps.txt') then
      Append(FileOutReplicate)
    Else
     begin
      ReWrite(FileOutReplicate);
      Write(FileOutReplicate,
            'RndSeed' + #9 +
            'InitCell' + #9);

      for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
        Write(FileOutReplicate,
              'CarCapParam' + IntToStr(i+1) + #9);

      Write(FileOutReplicate,
            'FitIdx' + #9 +
            'nLangs' + #9 +
            'r2 Rich' + #9 +
            'r2 RS' + #9 +
            'D');
     end;
   end;






  SetLength(SimVec, nThreads);
  SetLength(HandleVec, nThreads);
  For i:= 0 to nThreads-1 do
   begin
    SimVec[i]:= TSim.Create(True, True, True, True);
    SimVec[i].FreeOnTerminate:= True;
    SimVec[i].Priority:= tpLower;
    SimVec[i].SimId:= i;
    HandleVec[i]:= SimVec[i].Handle;

    SimVec[i].TotNumReps:= RepsPerThreads;

    SimVec[i].ModelParameters:= ModelParameters;

// Duplicate map data to each thread
    SimVec[i].MapShape:= MapShape;
    SimVec[i].MapShapeSmallRes:= MapShapeSmallRes;

    If TTextRec(FileOutReplicate).Mode <> 0 then
     begin
      SimVec[i].FileOutReplicate:= @FileOutReplicate;
      SimVec[i].FileOutReplicateWrite:= @FileOutReplicateWrite;
     end;

    SimVec[i].TotNumCells:= TotNumCells;
    SimVec[i].TotNumEnvVars:= TotNumEnvVars;
    SetLength(SimVec[i].MapCells, TotNumCells);
    SetLength(SimVec[i].CellSeq, TotNumCells);
    For c:= 0 to TotNumCells-1 do
     begin
      SimVec[i].MapCells[c].IdCell:= MapCells[c].IdCell;
      SimVec[i].MapCells[c].IdCell:= MapCells[c].IdCell;
      SimVec[i].MapCells[c].Coord.X:= MapCells[c].Coord.X;
      SimVec[i].MapCells[c].Coord.Y:= MapCells[c].Coord.Y;
      SimVec[i].MapCells[c].Area:= MapCells[c].Area;
      SimVec[i].CellSeq[c]:= CellSeq[c];

      SetLength(SimVec[i].MapCells[c].Environment, TotNumEnvVars);
      For d:= 0 to TotNumEnvVars-1 do
        SimVec[i].MapCells[c].Environment[d]:= MapCells[c].Environment[d];

      SetLength(SimVec[i].MapCells[c].Neighbors, Length(MapCells[c].Neighbors));
      For d:= 0 to Length(MapCells[c].Neighbors)-1 do
        SimVec[i].MapCells[c].Neighbors[d]:= @SimVec[i].MapCells[MapCells[c].Neighbors[d].IdCell];
     end;

    SimVec[i].SetUpModel(False);
   end;

  // If the total number of replicates is not a multiple of the number of threads
  j:= RepsPerThreads * nThreads;
  i:= 0;
  while j <> nReps do
   begin
    // Add additional reps to the first threads, until the total number of reps add up to nReps
    SimVec[i].TotNumReps:= SimVec[i].TotNumReps + 1;
    i:= i + 1;
    j:= j + 1;
   end;

  For i:= 0 to nThreads-1 do
    SimVec[i].Start;

  WaitForMultipleObjects(nThreads, Pointer(HandleVec), True, INFINITE);











  // Average the fit index calculated by each thread
  Result:= 0;
  s:= 0;

  AvgFitStats.nLangs:= 0;

  AvgFitStats.KS:= 0;
  AvgFitStats.RegResAvgRS.r2:= 0;
  AvgFitStats.RegResAvgRS.r2Adj:= 0;

  AvgFitStats.RegResAvgRich.r2:= 0;
  AvgFitStats.RegResAvgRich.r2Adj:= 0;

  For i:= 0 to Length(SimVec)-1 do
   begin
    AvgFitStats.nLangs:= AvgFitStats.nLangs + (SimVec[i].AvgFitStats.nLangs * SimVec[i].TotNumReps);

    AvgFitStats.KS:= AvgFitStats.KS + (SimVec[i].AvgFitStats.KS * SimVec[i].TotNumReps);
    AvgFitStats.RegResAvgRS.r2:= AvgFitStats.RegResAvgRS.r2 + (SimVec[i].AvgFitStats.RegResAvgRS.r2 * SimVec[i].TotNumReps);
    AvgFitStats.RegResAvgRS.r2Adj:= AvgFitStats.RegResAvgRS.r2Adj + (SimVec[i].AvgFitStats.RegResAvgRS.r2Adj * SimVec[i].TotNumReps);

    AvgFitStats.RegResAvgRich.r2:= AvgFitStats.RegResAvgRich.r2 + (SimVec[i].AvgFitStats.RegResAvgRich.r2 * SimVec[i].TotNumReps);
    AvgFitStats.RegResAvgRich.r2Adj:= AvgFitStats.RegResAvgRich.r2Adj + (SimVec[i].AvgFitStats.RegResAvgRich.r2Adj * SimVec[i].TotNumReps);

    Result:= Result + (SimVec[i].FitIndex * SimVec[i].TotNumReps);
    s:= s + SimVec[i].TotNumReps;

    //SimVec[i].Free;
   end;

  AvgFitStats.nLangs:= AvgFitStats.nLangs / s;

  AvgFitStats.KS:= AvgFitStats.KS / s;
  AvgFitStats.RegResAvgRS.r2:= AvgFitStats.RegResAvgRS.r2 / s;
  AvgFitStats.RegResAvgRS.r2Adj:= AvgFitStats.RegResAvgRS.r2Adj / s;

  AvgFitStats.RegResAvgRich.r2:= AvgFitStats.RegResAvgRich.r2 / s;
  AvgFitStats.RegResAvgRich.r2Adj:= AvgFitStats.RegResAvgRich.r2Adj / s;

  Result:= Result / s;
  FitIndex:= Result;






//
// These are the final results for each parameter combination, averaged among replicates of a thread and among threads
//  WriteLn('');
  WriteLn('');
  Write('  --->   ',
        ModelParameters.InitCell:0, '  ');

  for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
    Write(ModelParameters.CarryingCapParams[i]:0:3, '  ');

  Write(FitIndex:0:3, '  ',
        AvgFitStats.nLangs:0:0, '  ',
        AvgFitStats.RegResAvgRich.r2:0:3, '  ',
        AvgFitStats.RegResAvgRS.r2:0:3, '  ',
        AvgFitStats.KS:0:3);




  // These are the final results, averaged among replicates, for a given parameter combination
  AssignFile(FileOutSim, ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - Averages.txt');
  If FileExists(ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - Averages.txt') then
    Append(FileOutSim)
  Else
   begin
    ReWrite(FileOutSim);

    Write(FileOutSim, 'CarCapFunc' + #9);

    for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
      Write(FileOutSim,
            'CarCapParam' + IntToStr(i+1) + #9);

    Write(FileOutSim,
          'FitIdx' + #9 +
          'nLangs' + #9 +
          'r2 Rich' + #9 +
          'r2 RS' + #9 +
          'D');
   end;

  WriteLn(FileOutSim);

  if ModelParameters.ModelFunction = TMFPower then
    Write(FileOutSim, 'Power', #9) Else
  if ModelParameters.ModelFunction = TMFLogistic then
    Write(FileOutSim, 'Logistic', #9) Else
  if ModelParameters.ModelFunction = TMFExponential then
    Write(FileOutSim, 'Exp', #9);

  for i := 0 to Length(ModelParameters.CarryingCapParams)-1 do
    Write(FileOutSim,
          ModelParameters.CarryingCapParams[i]:0:7, #9);

  Write(FileOutSim,
        FitIndex:0:5, #9,
        AvgFitStats.nLangs:0:3, #9,
        AvgFitStats.RegResAvgRich.r2:0:5, #9,
        AvgFitStats.RegResAvgRS.r2:0:5, #9,
        AvgFitStats.KS:0:5);

  CloseFile(FileOutSim);
//


  // This should never happen!
  If IsNan(Result) then
   begin
    CloseFile(FileOutReplicate);
    WriteLn('Result is NaN');
    ReadLn(Result);
   end;






  SimVec:= Nil;

  If TTextRec(FileOutReplicate).Mode <> 0 then
   begin
    DeleteCriticalSection(FileOutReplicateWrite);
    CloseFile(FileOutReplicate);
   end;
 end;

Function TSim.Sample(List: TIntVector; n: Integer = -1; WithReplacement: Boolean = False): TIntVector;
 Function IntRND(Min, Max: TFloat): Integer;
  begin
   Result := Trunc(((Max - Min + 1) * (Random) ) + Min);
  end;
 var
  i: Integer;
  SeqOrder: TDblMatrix;
 begin
  If List = Nil then
   begin
    Result:= Nil;
    Exit;
   end;

  If n = -1 then
    n:= Length(List);

  SetLength(Result, n);
  If WithReplacement then
   begin
    For i:= 0 to n-1 do
     begin
      Result[i]:= List[IntRND(0,Length(List)-1)];
     end;
   end
  Else
   begin
    If n > Length(List) then
      raise Exception.Create('Impossible to sample without replacement more than the population size');
    SetLength(SeqOrder, Length(List), 2);
    For i:= 0 to Length(SeqOrder)-1 do
     begin
      SeqOrder[i,0]:= i;
      SeqOrder[i,1]:= Random;
     end;
    SeqOrder:= QuickSort(SeqOrder, 1);

    SetLength(Result, n);
    For i:= 0 to n-1 do
     begin
      Result[i]:= List[Trunc(SeqOrder[i,0])];
     end;
   end;
 end;

Function TSimParallel.ExecToGibbsSamp(ProposalParams: TDblVector): Double;
 var
  i: Integer;
 begin
  for i:= 0 to Length(ProposalParams)-1 do
    ModelParameters.CarryingCapParams[i]:= ProposalParams[i];

  Execute(ModelParameters,
          tnReps,
          tnThreads);

  Result:= FitIndex;
 end;








































//MCMC search only used in GEB paper of 2017
Procedure TSimParallel.SearchGibbs(ModelParameters: TModelParameters;
                                   DisturbFuncs: TDisturbFuncs;
                                   var MCMCChain: TChain;
                                   nReps: Integer = 120;
                                   nThreads: Integer = -1);
 var
  StartValues: TDblVector;
  Priors: TPriors;
  i: Integer;
 begin
  SetLength(DisturbFuncs, Length(ModelParameters.CarryingCapParams));

  if ModelParameters.ModelFunction = TMFPower then
   begin
    // Carrying Capacity Scale Parameter
    DisturbFuncs[0].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[0].SampFuncPars, 2);
    DisturbFuncs[0].SampFuncPars[0]:= 0;
    DisturbFuncs[0].SampFuncPars[1]:= 0.05;

    // Carrying Capacity Slope-Rate Parameter
    DisturbFuncs[1].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[1].SampFuncPars, 2);
    DisturbFuncs[1].SampFuncPars[0]:= 0;
    DisturbFuncs[1].SampFuncPars[1]:= 0.05;
   end Else

  if ModelParameters.ModelFunction = TMFLogistic then
   begin
    // x0 is the value of precipitation in which the curve is in its midpoint (point of inflection of the curve)
    DisturbFuncs[0].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[0].SampFuncPars, 2);
    DisturbFuncs[0].SampFuncPars[0]:= 0;
    DisturbFuncs[0].SampFuncPars[1]:= 10;

    // s is the steepness of the curve ; because s is a very small number, we estimate log10(s)
    DisturbFuncs[1].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[1].SampFuncPars, 2);
    DisturbFuncs[1].SampFuncPars[0]:= 0;
    DisturbFuncs[1].SampFuncPars[1]:= 0.01;

    // m is the maximum carrying capacity (saturation level), in units of carrying capacity
    DisturbFuncs[2].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[2].SampFuncPars, 2);
    DisturbFuncs[2].SampFuncPars[0]:= 0;
    DisturbFuncs[2].SampFuncPars[1]:= 1;
   end Else

  if ModelParameters.ModelFunction = TMFExponential then
   begin
    // a is a scaler
    DisturbFuncs[0].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[0].SampFuncPars, 2);
    DisturbFuncs[0].SampFuncPars[0]:= 0;
    DisturbFuncs[0].SampFuncPars[1]:= 0.01;

    // b is a multiplier of precipitation, which is in the exponent of euler number (natural exponent)
    DisturbFuncs[1].RndDistrib:= TRndNorm;
    SetLength(DisturbFuncs[1].SampFuncPars, 2);
    DisturbFuncs[1].SampFuncPars[0]:= 0;
    DisturbFuncs[1].SampFuncPars[1]:= 0.005;
   end;

  SetLength(Priors, Length(ModelParameters.CarryingCapParams));

  Self.ModelParameters:= ModelParameters;
  tnReps:= nReps;
  tnThreads:= nThreads;

  SetLength(StartValues, Length(ModelParameters.CarryingCapParams));
  for i:= 0 to Length(ModelParameters.CarryingCapParams)-1 do
    StartValues[i]:= ModelParameters.CarryingCapParams[i];

  GibbsSamp:= TGibbs.Create(StartValues,             // Initial values for the MCMC
                            ExecToGibbsSamp,
                            Priors,
                            DisturbFuncs,
                            ModelParameters.OutPutFolder + '!Sim' + ModelParameters.SimName + ' - MCMC Chain.txt',
                            100000,
                            0);
 end;

end.
