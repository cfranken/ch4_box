function [ obs ] =replaceMCF(obs,St)

sNames = {'alt','brw','kum','mlo','smo','cgo','psa','spo','mhd','thd','ush','sum'};
sLat   = [ 82.5, 71.3, 19.5, 19.5, 14.3,-40.7,-64.6,  -90,  53,   41, -54.8, 72.6];
sCol   = [    9,   13,    21,   23,   25,   27,   31,   33,  15,   17,    29,   11];
noaa = load('data/obs/mcf/NOAA/combined/HATS_global_MC_mat.txt');

for i=1:length(sCol)
    a = isfinite(noaa(:,sCol(i)));
    lRecord(i) = sum(a);
end

NH = find(sLat>10 & lRecord>200);
SH = find(sLat<-10 & lRecord>200);


tt = datenum(noaa(:,1), noaa(:,2), 15);
TD = (St(2)-St(1))/2.

for i=1:length(St)
    wo = find(abs(tt-St(i))<TD);
    disp(i)
    disp(length(wo))
    if (length(wo)>11 & St(i)>datenum(1991,12,1))
        obs.nh_mcf(i) = nanmean(nanmean(noaa(wo,sCol(NH))));
        obs.sh_mcf(i) = nanmean(nanmean(noaa(wo,sCol(SH))));
        %if nanstd(nanmean(noaa(wo,sCol(NH)))) > obs.nh_mcf_err(i)
            obs.nh_mcf_err(i) = max(nanstd(nanmean(noaa(wo,sCol(NH)))),0.5/100*obs.nh_mcf(i));
        %end
        %if nanstd(nanmean(noaa(wo,sCol(SH)))) > obs.sh_mcf_err(i)
        
            obs.sh_mcf_err(i) = max(nanstd(nanmean(noaa(wo,sCol(SH)))),0.5/100*obs.sh_mcf(i));
        %end
        %obs.nh_mcf_err(i) = nanstd(nanmean(noaa(wo,sCol(NH))));
        %obs.sh_mcf_err(i) = nanstd(nanmean(noaa(wo,sCol(SH))));
        %if obs.sh_mcf_err(i)==0
        %    obs.sh_mcf_err(i)= 0.1;
        %end
    else
        obs.nh_mcf(i) = NaN;
        obs.sh_mcf(i) = NaN;
        obs.sh_mcf_err(i) = 200;
        obs.nh_mcf_err(i) = 200;
    end
    
end
