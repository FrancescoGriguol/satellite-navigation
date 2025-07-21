function [SAT,SatID] = brdc_get_sc_data(brdc,SatID,brdcEphem,constID)

i_el = [];

for i = 1:length(SatID)
    if nnz(brdc.(constID).SatelliteID == SatID(i)) == 0
        disp(['SatID =', num2str(SatID(i)),' not found in RINEX file'])
        i_el = [i_el,i];
        continue
    end
    str = string(['SatID_', num2str(SatID(i),'%i')]);
    SAT.(str) = brdcEphem(brdcEphem.SatelliteID == SatID(i),:);

end

SatID(i_el) = [];

end