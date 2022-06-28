%
%%

%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Paths ==========================

OutputDir = [BaseDir '\data\preprocess\yeast_tc\E1R1\Exp1_rep1_0min_im1\STL1-TMR'];

% ========================== Metadata ==========================

rna_targ = 'STL1';
rna_probes = 'TMR';
rna_moltype = 'mRNA';
cell_type = 'BY4741';
cell_species = 'Saccharomyces cerevisiae';

% ========================== Test ==========================

addpath('./core');

Main_RNASpotsMetaAdd(OutputDir, '-probetype', rna_probes, '-target', rna_targ, '-targettype', rna_moltype, '-species', cell_species, '-celltype', cell_type);