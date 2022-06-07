%

%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Paths ==========================

InputDir = [BaseDir '\img\yeast_full\Exp1_rep1'];
CellsegDir = [BaseDir '\data\cell_seg\yeast_full\E1R1'];
OutputDir = [BaseDir '\data\preprocess\yeast_tc\E1R1'];

% ========================== Image Info ==========================

ch_count = 4;
light_ch = 4;
rna_ch = [1 2];
rna_targ = {'STL1' 'CTT1'};
rna_probes = {'TMR' 'CY5'};
rna_moltype = {'mRNA' 'mRNA'};
cell_type = 'BY4741';
cell_species = 'Saccharomyces cerevisiae';
z_trim = 3;

tmax = [400 400];
ttune_mad = [0.0 0.0];
ttune_winsize = [15 15];

% ========================== Get a list of all images ==========================

DirContents = dir(InputDir);
fcount = size(DirContents,1);

j = 1;
FileList = cell(fcount, 1);
for i = 1:fcount
    fname = DirContents(i,1).name;
    if ~DirContents(i,1).isdir && endsWith(fname, ".tif")
        FileList{j,1} = fname;
        j = j+1;
    end
end
fcount = j-1;

% ========================== Run ==========================

for i = 1:fcount
    fname = FileList{i,1};
    fprintf("Now processing %s (%d of %d)\n", fname, i, fcount);
    
    tifpath = [InputDir filesep fname];
    [~, imgname, ~] = fileparts(tifpath);
    cellseg_path = [CellsegDir filesep 'Lab_' imgname '.mat'];
    if ~isfile(cellseg_path)
        cellseg_path = [];
    end
    
    outdir = [OutputDir filesep imgname];
    probe_ch_count = size(rna_ch,2);
    for j = 1:probe_ch_count
        c_suffix = [rna_targ{1,j} '-' rna_probes{1,j}];
        choutdir = [outdir filesep c_suffix];
        chname = [imgname '_' c_suffix];
        
        %RN doesn't take control info...
        spotsrun = Main_RNASpots(chname, tifpath, rna_ch(1,j), light_ch, ch_count,...
        choutdir, 10, tmax(1,j), z_trim, cellseg_path, [], 0, 0, ttune_winsize(1,j), ttune_mad(1,j), false);
    
        spotsrun.type_probe = rna_probes{1,j};
        spotsrun.type_target = rna_targ{1,j};
        spotsrun.type_species = cell_species;
        spotsrun.type_cell = cell_type;
        spotsrun.type_targetmol = rna_moltype{1,j};
        spotsrun.saveMe();
    end
    
end

