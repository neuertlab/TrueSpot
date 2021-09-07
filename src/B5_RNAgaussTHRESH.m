function [PARTMR] = B5_RNAgaussTHRESH(TMR3Dorg,Lab,PAR,mm,Nuc)
clear PARTMR  
PARTMR = PAR;
nn = size(PARTMR.xfit,2);
PARTMR.nucRNA = NaN(mm,nn);
PARTMR.cytoRNA = NaN(mm,nn);
A = size(TMR3Dorg,1);
zz = size(TMR3Dorg,3);  %Image size
ps = 65;
for i = 1:mm
    i
    %Find each nucleus
    k1 = Lab == i ; %== q;%sieve out dots in cell q.
    k1=uint16(k1);
    k3=regionprops(k1,'BoundingBox','Area');
    k4 = k3.BoundingBox; %create the rectangular box around the cell. 
    X0=round(k4(1))-4;
    Y0=round(k4(2))-4;
    X1=round(k4(1)+ k4(3))+4;
    Y1=round(k4(2)+ k4(4))+4;
    if X0 < 1;
        X0 = 1;
    else;
    end;
    if Y0 < 1;
        Y0 = 1;
    else;

    end;
    if X1 > A;
        X1 = A;
    else;
    end;
    if Y1 > A;
        Y1 = A;
    else;
    end;
    k2 = k1(Y0:Y1,X0:X1);

    %% Separate Cell, Nuclear and cyto plasmic RNA
    clear smNuc

    for q=1:zz;
         smNuc(:,:,q) = immultiply(double(Nuc(Y0:Y1,X0:X1,q)),double(k2)); % Nucleus object
    end;
    %for all RNA spots
    n = 1;
    for n = 1:nn
        %n
      if ~isnan(PARTMR.xfit(i,n))
            PARTMR.nucRNA(i,n) = smNuc(round(PARTMR.xfit(i,n)),round(PARTMR.yfit(i,n)),PARTMR.zinit(i,n)) > 0; % Nuclear RNA
            if PARTMR.nucRNA(i,n) ~= 0;
                for v=1:(A./2)
                    if smNuc(round(PARTMR.xfit(i,n)+v),round(PARTMR.yfit(i,n)),PARTMR.zinit(i,n)) == 0
                        xloc = v;
                        break;
                    elseif smNuc(round(PARTMR.xfit(i,n)-v),round(PARTMR.yfit(i,n)),PARTMR.zinit(i,n)) == 0
                        xloc = v;
                        break;
                    end
                end
                for v =1:(A./2)
                    if smNuc(round(PARTMR.xfit(i,n)),round(PARTMR.yfit(i,n)+v),PARTMR.zinit(i,n)) == 0
                        yloc = v;
                        break;
                    elseif smNuc(round(PARTMR.xfit(i,n)),round(PARTMR.yfit(i,n)-v),PARTMR.zinit(i,n)) == 0
                        yloc = v;
                        break;
                    end
                end
                if xloc < yloc
                    PAR.distRNA(i,n) = xloc*ps;  %the result of the distance in nm of all nuclear RNA
                elseif yloc < xloc
                    PAR.distRNA(i,n) = yloc*ps;
                end 
                PARTMR.szNUC(i,n) = nnz(smNuc(:,:,PARTMR.zinit(i,n)))*ps.^2;
                PARTMR.normdist(i,n) = PARTMR.distRNA(i,n)./sqrt(PARTMR.szNUC(i,n)./pi());
            end
        PARTMR.cytoRNA(i,n) = smNuc(round(PARTMR.xfit(i,n)),round(PARTMR.yfit(i,n)),PARTMR.zinit(i,n)) == 0; % Cytoplasmic
        %Now find the cytoplasmic ellipticity
        clear len1 len2
        len1 = size(smNuc(:,:,PARTMR.zinit(i,n)),1)
        len2 = round(len1./2)
        if len1 > (size(smNuc,1));
            len1 = (size(smNuc,1))
        else
        end;
        if len2 > (size(smNuc,2));
            len2 = (size(smNuc,2))
        else
        end;
        
        for v = 1:(A./2)
            %v
            if smNuc(len1+v-1,len2,PARTMR.zinit(i,n)) == 0
            mxloc = abs(v); 
            break;
            end  
        end

        for v = 1:(A./2)
            
            if smNuc(len1-v+1,len2,PARTMR.zinit(i,n)) == 0
                pxloc = abs(v);
                break;
            end
        end

        for v = 1:(A./2)
         
            if smNuc(len1,len2+v-1,PARTMR.zinit(i,n)) == 0
                myloc = abs(v);
                break;
            end
        end

        for v = 1:(A./2)
            
            if smNuc(len1,len2-v+1,PARTMR.zinit(i,n)) == 0
                pyloc = abs(i);
                break
            end
        end

        PARTMR.nucR(i,n) = (mxloc+pxloc)./(myloc+pyloc);        
        end
    end;

    end
end
