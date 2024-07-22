%
%%
classdef RNACloud
    
    properties
        cloud_box; %Coordinates of cloud, relative to its cell
        cloud_mask;
        cloud_data; 
        
        is_nuc;
        is_cyto;
        nascent_flag = false;
        
        cloud_volume; %Just the number of pix/vox cloud encompasses.
        max_intensity; %TODO Add
        total_intensity;
    end
    
    methods
        
        %%
        function obj = recalculateVolume(obj)
            if isempty(obj.cloud_mask)
                obj.cloud_volume = 0.0;
                return; 
            end
            obj.cloud_volume = sum(obj.cloud_mask, 'all');
        end
        
    end
    
    methods (Static)
        
        %%
        function rna_cloud = newRNACloud()
            rna_cloud = RNACloud;
            rna_cloud.cloud_box = SingleCell.generateRecPrismStruct(1,1,1,1,1,1);
            rna_cloud.cloud_mask = [];
            rna_cloud.is_nuc = false;
            rna_cloud.is_cyto = false;
            rna_cloud.cloud_volume = 0.0;
            rna_cloud.total_intensity = 0.0;
            rna_cloud.cloud_data = [];
        end
        
    end
    
end