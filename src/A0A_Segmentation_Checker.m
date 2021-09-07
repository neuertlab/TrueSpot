function A0A_Segmentation_Checker(dates0,times0,max_tps,max_imnum,CY5_name,AF594_name,TMR_name,nucTH)

for i = 1:size(dates0,2)
    for j = 1:max_tps
        for k = 1:max_imnum
            try
                filename = ['Lab_' dates0{i} '_' times0{i}{j} '_mESC_' CY5_name{i}(2:size(CY5_name{i},2)-2) AF594_name{i}(2:size(AF594_name{i},2)-2) TMR_name{i}(2:size(TMR_name{i},2)-3) '_img' num2str(k)]
                load(filename,'cells','trans_plane');
                figure(101); clf; title([dates0{i} ' ' times0{i}{j} ' img ' num2str(k)])
                subplot(1,2,1); imshow(cells,[]); %title([dates0{i} ' ' times0{i}{j} ' img ' num2str(k)])
                subplot(1,2,2); imshow(trans_plane,[]); %title([dates0{i} ' ' times0{i}{j} ' img ' num2str(k)])
                pause(3)
            catch
                
            end
        end
    end
end
