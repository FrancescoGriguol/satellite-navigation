function sp3 = sp3_compute_sc_pos_eci(sp3,Greg_time, constID)

no_sat = height(sp3.data(:,1))/length(Greg_time); % number of satellites
sp3.(strcat(constID,'_ECI')) = sp3.data;   % preallocation
sp3.(strcat(constID,'_ECI')){:,5:7} = NaN(height(sp3.data(:,1)),3);

for i = 1:length(Greg_time)
    for j = 1 : no_sat
        r_ecef = table2array(sp3.data(j+no_sat*(i-1),5:7));
        r_eci = ecef2eci(Greg_time(i,:),r_ecef);
        sp3.(strcat(constID,'_ECI')){j+no_sat*(i-1),5:7} = r_eci';
    end
end

end