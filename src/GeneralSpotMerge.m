%%% Below removes duplicates of the same
%%% transcript. Goes through every transcript,
%%% finds distance between it and other
%%% transcripts, and eliminates others

for p = 1:size(PARcy5_mid.xinit,2)
    dists = zeros(1,size(size(PARcy5_mid.xinit,2),3)); %stores distance from transcript to others
    if not(isnan(PARcy5_mid.xinit(m,p)))
        for pp = 1:size(PARcy5_mid.xinit,2)    %find distance squared
            dists(1,pp,1) = (PARcy5_mid.xinit(m,p)- PARcy5_mid.xinit(m,pp))^2;
            dists(1,pp,2) = (PARcy5_mid.yinit(m,p)- PARcy5_mid.yinit(m,pp))^2;
            dists(1,pp,3) = ((PARcy5_mid.zinit(m,p)- PARcy5_mid.zinit(m,pp))*z_adjust)^2;
        end
        distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
        for pp = 1:size(PARcy5_mid.xinit,2)
            if p ~= pp & distF(pp) < same_dist
                if PARcy5_mid.TotExpInt(m,pp) < PARcy5_mid.TotExpInt(m,p) %Only keep the brighter spot
                    PARcy5_mid.xinit(m,pp) = NaN;
                    PARcy5_mid.yinit(m,pp) = NaN;
                    PARcy5_mid.zinit(m,pp) = NaN;
                    PARcy5_mid.TotExpInt(m,p) = PARcy5_mid.TotExpInt(m,p)+PARcy5_mid.TotExpInt(m,pp);
                    PARcy5_mid.TotExpInt(m,pp) = NaN;
                else
                    PARcy5_mid.xinit(m,p) = NaN;
                    PARcy5_mid.yinit(m,p) = NaN;
                    PARcy5_mid.zinit(m,p) = NaN;
                    PARcy5_mid.TotExpInt(m,pp) = PARcy5_mid.TotExpInt(m,p)+PARcy5_mid.TotExpInt(m,pp);
                    PARcy5_mid.TotExpInt(m,p) = NaN;
                    break
                end
            end
        end
    end
end
