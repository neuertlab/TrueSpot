%
%%
% ========================== Params ==========================
addpath('./core');

ImgBaseDir = 'C:\Users\hospelb\labdata\imgproc\img';

ScanDir = [ImgBaseDir '\histones_febc'];
do_recursive = false;

% ========================== Go ==========================

dofordir(ScanDir, do_recursive);

% ========================== Scan Func ==========================

function dofordir(dirpath, r)
    state_obj = struct('is_okay', true);
    dircontents = dir(dirpath);
    fcount = size(dircontents,1);
    for i = 1:fcount
        file_rec = dircontents(i);
        if file_rec.isdir
            if r
                if strcmp(file_rec.name, '.')
                    continue;
                elseif strcmp(file_rec.name, '..')
                    continue;
                else
                    dofordir([file_rec.folder filesep file_rec.name], r);
                end
            end
        else
            filepath = [file_rec.folder filesep file_rec.name];
            doathing(filepath, state_obj);
        end
        if ~state_obj.is_okay; return; end
    end
end

% ========================== Do Func ==========================

function doathing(filepath, state_obj)
    if ~endsWith(filepath, '.tif'); return; end
    [~, fname, ~] = fileparts(filepath);
    iname_base = ['histonesc_' fname];
    
    for i = 2:4
        target = '';
        if i == 2
            target = 'Xist'; 
            probe = 'CY5';
        elseif i == 3
            target = 'Tsix';
            probe = 'TMR';
        elseif i == 4
            target = 'Histone'; 
            probe = 'AF488';
        end
        
        fprintf("%s\t", [iname_base '_' target]);
        fprintf("%s\t", ['/img/histones_febc/' fname '.tif']);
        fprintf("%d\t5\t1\t5\t", i);
        
        if i == 4
            if contains(fname, 'H3K36me3'); target = 'H3K36me3'; end
            if contains(fname, 'H3K4me2'); target = 'H3K4me2'; end
        end
        fprintf("%s\t", ['/data/preprocess/histones_febc/' fname '/' target '-' probe '/' fname '-' target '-' probe '_all_3d']);
        fprintf("mESC\t%s\t%s\t", target, probe);
        if i == 4
            fprintf("Histone Mark\t");
        else
            fprintf("lncRNA\t");
        end
        fprintf("Mus musculus\t");
        fprintf("%s\t", ['/data/bigfish/histones_febc/' fname '/' target '-' probe '/BIGFISH_' fname '-' target '-' probe]);
        fprintf("3\t.\t/data/cell_seg/histones_febc\t%s\t", fname);
        fprintf("65\t65\t300\t");
        if i == 4
            fprintf("160\t160\t300\t");
        else
            fprintf("195\t195\t310\t");
        end
        fprintf("0\t0\t.");
        fprintf("\n");
    end
    
    
end

% function doathing(filepath, state_obj)
%     if ~endsWith(filepath, '.tif'); return; end
%     if endsWith(filepath, '_C0.tif'); return; end
%     if endsWith(filepath, '_C1.tif'); return; end
%     if endsWith(filepath, '_C2.tif'); return; end
%     if endsWith(filepath, '_C3.tif'); return; end
%     
%     [parentdir, iname, ~] = fileparts(filepath);
%     [grandparentdir, parentdir, ~] = fileparts(parentdir);
%     parentdir = replace(parentdir, filesep, '/');
%     
%     if endsWith(iname, '_rnasum')
%         fprintf("%s\t", iname);
%         [~, parentdir, ~] = fileparts(grandparentdir);
%         parentdir = replace(parentdir, filesep, '/');
%         fprintf("%s\t1\t0\t0\t1\t", ['/img/munsky_lab/' parentdir '/rnasum/' iname '.tif']);
%         fprintf("%s\tHeLa\tmRNA\tGFPwCY5\tmRNA\tHomo sapiens\t", ['/data/preprocess/munsky_lab/' parentdir '/' iname '/mRNA-GFPwCY5/' iname '-mRNA-GFPwCY5_all_3d']);
%         fprintf("%s\t", ['/data/bigfish/munsky_lab/' parentdir '/' iname '/mRNA-GFPwCY5/BIGFISH_' iname '-mRNA-GFPwCY5']);
%     else
%         fprintf("%s_GFP\t", iname);
%         fprintf("%s\t2\t0\t1\t4\t", ['/img/munsky_lab/' parentdir '/' iname '.tif']);
%         fprintf("%s\tHeLa\tmRNA\tGFP\tmRNA\tHomo sapiens\t", ['/data/preprocess/munsky_lab/' parentdir '/' iname '/mRNA-GFP/' iname '-mRNA-GFP_all_3d']);
%         fprintf("%s\t", ['/data/bigfish/munsky_lab/' parentdir '/' iname '/mRNA-GFP/BIGFISH_' iname '-mRNA-GFP']);
%         
%         fprintf("3\t.\t.\t.\t");
%         fprintf("160\t160\t500\t120\t120\t350\t");
%         fprintf("0\t0\t.\n");
%         
%         fprintf("%s_CY5\t", iname);
%         fprintf("%s\t4\t0\t1\t4\t", ['/img/munsky_lab/' parentdir '/' iname '.tif']);
%         fprintf("%s\tHeLa\tmRNA\tCY5\tmRNA\tHomo sapiens\t", ['/data/preprocess/munsky_lab/' parentdir '/' iname '/mRNA-CY5/' iname '-mRNA-CY5_all_3d']);
%         fprintf("%s\t", ['/data/bigfish/munsky_lab/' parentdir '/' iname '/mRNA-CY5/BIGFISH_' iname '-mRNA-CY5']);
%     end
%     fprintf("3\t.\t.\t.\t");
%     fprintf("160\t160\t500\t120\t120\t350\t");
%     fprintf("0\t0\t.\n");
%     
% end

% function doathing(filepath, state_obj)
%     if ~endsWith(filepath, '.tif'); return; end
%     if endsWith(filepath, '_C0.tif'); return; end
%     if endsWith(filepath, '_C1.tif'); return; end
%     if endsWith(filepath, '_C2.tif'); return; end
%     if endsWith(filepath, '_C3.tif'); return; end
%     if endsWith(filepath, '_rnasum.tif'); return; end
%     
%     fprintf("Found compatible tif: %s\n", filepath);
%     [channels, ~] = LoadTif(filepath, 4, [2 4], 1);
%     ch_gfp = channels{2,1};
%     ch_cy5 = channels{4,1};
%     ch_sum = ch_gfp + ch_cy5;
%     
%     highest_intensity = max(ch_sum, [], 'all');
%     if highest_intensity <= 65536
%         ch_sum = uint16(ch_sum);
%     else
%         fprintf("Sum would clip 16 bits. Saving as 32 bit...\n");
%         ch_sum = uint32(ch_sum);
%     end
%     
%     [dirpath, fname, ~] = fileparts(filepath);
%     outpath = [dirpath filesep fname '_rnasum.tif'];
%     if isfile(outpath) 
%         delete(outpath);
%     end
%     saveastiff(ch_sum, outpath);
%     
% end
