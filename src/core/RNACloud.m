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
        
        %%
        function pkg = packageForSave(obj)
            pkg = struct();
            pkg.cloud_box = obj.cloud_box;
            pkg.cloud_mask = obj.cloud_mask;
            pkg.cloud_data = obj.cloud_data;
            pkg.is_nuc = obj.is_nuc;
            pkg.is_cyto = obj.is_cyto;
            pkg.nascent_flag = obj.nascent_flag;
            pkg.cloud_volume = obj.cloud_volume;
            pkg.max_intensity = obj.max_intensity;
            pkg.total_intensity = obj.total_intensity;
        end
    
    end
    
    methods (Static)

        %%
        function rna_cloud = readFromSavePackage(pkg)
            rna_cloud = RNACloud;
            if isempty(pkg); return; end
            rna_cloud.cloud_box = pkg.cloud_box;
            rna_cloud.cloud_mask = pkg.cloud_mask;
            rna_cloud.cloud_data = pkg.cloud_data;
            rna_cloud.is_nuc = pkg.is_nuc;
            rna_cloud.is_cyto = pkg.is_cyto;
            rna_cloud.nascent_flag = pkg.nascent_flag;
            rna_cloud.cloud_volume = pkg.cloud_volume;
            rna_cloud.max_intensity = pkg.max_intensity;
            rna_cloud.total_intensity = pkg.total_intensity;
        end

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