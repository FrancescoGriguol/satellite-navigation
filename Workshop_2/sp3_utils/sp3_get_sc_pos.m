function [SAT,SatID] = sp3_get_sc_pos(sp3,SatID)
i_el = [];
% SP3
for i = 1 : length(SatID)
    if nnz(sp3.data.PRN == SatID(i))==0 
        %se ==0 allora significa che nella colonna di nostro interesse c'è
        % almeno uno dei datio che ci interessa (0=no, 1=si)
        disp(['SatID = ', num2str(SatID(i)), ' not found in SP3'])
        i_el = [i_el, i];
        continue
    end
    str = string(['SatID_', num2str(SatID(i),'%i')]);
    idx_r = sp3.data.PRN == SatID(i);
    GS = sp3.data( idx_r,sp3.col.X:sp3.col.Z ) ;
    SAT.(str) = table2array(GS);
end
SatID(i_el) = [];
end