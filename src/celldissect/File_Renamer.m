[images1,folders1] = Get_Dir_All('C:\Users\keslerbk\Dropbox (VU Basic Sciences)\CellSegmentationCodes\Yeast Analysis\Images .4M',1)
images1{1}
for i = 1:size(images1,2)
if strfind(images1{i},'_img_MMImages.ome.tif')
    start_mod = strfind(images1{i},'_img_MMImages.ome.tif');
    num_ind = strfind(images1{i},'\');
    num_ind = num_ind(size(num_ind,2))-1;
    img_num = images1{i}(num_ind);
    new_string = images1{i};
    new_string(start_mod:size(images1{i},2)+2) = ['_img_' img_num '_MMImages.ome.tif'];
    movefile(images1{i},new_string);
end
end
    
[images1,folders1] = Get_Dir_All('C:\Users\keslerbk\Desktop\Microscopy\20191030',1)
images1{1}
for i = 1:size(images1,2)
new_string = insertAfter(images1{i},'F1-2-1_Tsix','5')
    movefile(images1{i},new_string);
end
