function [TSiteN_Tsix,TSitePos_Tsix] = A1C_Transc_Site(All_spots_Tsix,Tsix_thres)

TSiteN_Tsix = 0;
TSitePos_Tsix = NaN(1,2,3);
%Tsix_thres = 3*median_int_Tsix;

for i = 1:size(All_spots,1)             %determine how many transcription sites there are
    TSite_temp = 0;
    TSite_temp(1,1:sum(All_spots_Tsix(i,:) > Tsix_thres),1:3) =  All_positions_Tsix(i,find(All_spots_Tsix(i,:) > Tsix_thres),:);
    sites = 1:sum(TSite_temp(1,:,1)>0);
    counter1 = 1;       %for looking at each site
    while counter1 <= size(sites,2) 
        diffsSq = 0; %column is specififc transc site compared to, 3rd d is x,y,z,
        dists = 0;
        for k = 1:size(sites,2)
            for m = 1:3
                diffsSq(1,k,m) = (TSite_temp(1,sites(k),m)- TSite_temp(1,sites(counter1),m))^2; %Find the Square difference
            end
            dists(1,k) = sqrt(diffsSq(1,k,1)+diffsSq(1,k,2)+diffsSq(1,k,3));
        end
        mults = find(dists < 10);    %Where the multiples are
        mults(mults<=counter1) = [];
        sites(mults) = [];
        counter1 = counter1+1;
        dists;
    end
    TSite_temp = TSite_temp(1,sites,:);
    TSitePos_Tsix(i,1:size(TSite_temp,2),1:3) = TSite_temp;
    TSiteN_Tsix(i,1) = sum(TSite_temp(:,:,1)>0);
end
max(TSiteN_Tsix)