%%
set(0,'defaultAxesFontSize',25)
clear; close all;
addpath(genpath(cd))

%%
global Seqdir TAXAdir OTUdir;

Seqdir='SEQanalysis/excludeNatural/M_OTUtableGreenGenes.csv';
TAXAdir='SEQanalysis/excludeNatural/M_TAXAtableGreenGenes.csv';
OTUdir='SEQanalysis/excludeNatural/M_UniqueSequences.csv';

PostProcessingSequences_synthetic
MetadataGenerator
RecipeGenerator
ExperimentalDataProcessing
%%
clear; 
ProcessedSeqdir='Postprocessed/processed_Sequences_synthetic.xlsx';
METAdir='Postprocessed/Metadata.xlsx';
CoalRecipedir='Postprocessed/CoalescenceRecipe.xlsx';
ProcessedSimdir='SimulationResult/result_rand/processed_simulations.xlsx';

global ProcessedSeq Metadata CoalRecipe CoalRecipe_natural ProcessedSim;
ProcessedSeq=readtable(ProcessedSeqdir);
Metadata=readtable(METAdir);
CoalRecipe=readtable(CoalRecipedir);
CoalRecipe_natural=readtable(CoalRecipedir,'Sheet', 2);
ProcessedSim=readtable(ProcessedSimdir)

AnalyzeCoalescence_synthetic
AnalyzeCoalescence_synthetic_simulations
AnalyzeCommunities_synthetic
%AnalyzeCommunities_synthetic_simulations



%%
ProcessedComdir='Analyzed/processed_Communities_synthetic.xlsx';
ProcessedCoaldir='Analyzed/processed_CoalescenceEvent_synthetic.xlsx';

Community_data=readtable(ProcessedComdir);
Metadata=readtable(METAdir);
Coalescence_data=readtable(ProcessedCoaldir);
