%%
set(0,'defaultAxesFontSize',25)
clear; close all;
addpath(genpath(cd))

%%
global Seqdir TAXAdir OTUdir;

Seqdir='SEQanalysis/onlyNatural/M_OTUtableGreenGenes.csv';
TAXAdir='SEQanalysis/onlyNatural/M_TAXAtableGreenGenes.csv';
OTUdir='SEQanalysis/onlyNatural/M_UniqueSequences.csv';

PostProcessingSequences_natural
%MetadataGenerator
%RecipeGenerator
%ExperimentalDataProcessing
%%
clear; 
ProcessedSeqdir='Postprocessed/processed_Sequences_natural.xlsx';
METAdir='Postprocessed/Metadata.xlsx';
CoalRecipedir='Postprocessed/CoalescenceRecipe.xlsx';
ProcessedSimdir='SimulationResult/result_rand/processed_simulations.xlsx';

global ProcessedSeq Metadata CoalRecipe CoalRecipe_natural ProcessedSim;
ProcessedSeq=readtable(ProcessedSeqdir);
Metadata=readtable(METAdir);
CoalRecipe=readtable(CoalRecipedir);
CoalRecipe_natural=readtable(CoalRecipedir,'Sheet', 2);
ProcessedSim=readtable(ProcessedSimdir)

AnalyzeCoalescence_natural
AnalyzeCommunities_natural


