
aaa = load('PLANAR_1coil_1port_sparam_scuff');
Nport = sqrt((size(aaa,2)-1)/2);

ff_scuff = squeeze(aaa(:,1))*1e9;
SP_scuff = zeros(Nport,Nport,length(ff_scuff));

for ii = 1:Nport
    rowshift = 2*Nport*(ii-1)+1;
    for jj = 1:Nport
        shift = rowshift + 2*(jj-1);
        SP_scuff(ii,jj,:) = aaa(:,shift+1) + 1j*aaa(:,shift+2);
    end
end
        
aaa = load('PLANAR_1coil_1port_zparam_scuff');
ZP_scuff = zeros(Nport,Nport,length(ff_scuff));

for ii = 1:Nport
    rowshift = 2*Nport*(ii-1)+1;
    for jj = 1:Nport
        shift = rowshift + 2*(jj-1);
        ZP_scuff(ii,jj,:) = aaa(:,shift+1) + 1j*aaa(:,shift+2);
    end
end

save('PLANAR_1coil_1port_scuff.mat', 'ff_scuff', 'SP_scuff', 'ZP_scuff');


aaa = load('PLANAR_1coil_2port_sparam_scuff');
Nport = sqrt((size(aaa,2)-1)/2);

ff_scuff = squeeze(aaa(:,1))*1e9;
SP_scuff = zeros(Nport,Nport,length(ff_scuff));

for ii = 1:Nport
    rowshift = 2*Nport*(ii-1)+1;
    for jj = 1:Nport
        shift = rowshift + 2*(jj-1);
        SP_scuff(ii,jj,:) = aaa(:,shift+1) + 1j*aaa(:,shift+2);
    end
end
        
aaa = load('PLANAR_1coil_2port_zparam_scuff');
ZP_scuff = zeros(Nport,Nport,length(ff_scuff));

for ii = 1:Nport
    rowshift = 2*Nport*(ii-1)+1;
    for jj = 1:Nport
        shift = rowshift + 2*(jj-1);
        ZP_scuff(ii,jj,:) = aaa(:,shift+1) + 1j*aaa(:,shift+2);
    end
end

save('PLANAR_1coil_2port_scuff.mat', 'ff_scuff', 'SP_scuff', 'ZP_scuff');



aaa = load('PLANAR_1coil_4port_sparam_scuff');
Nport = sqrt((size(aaa,2)-1)/2);

ff_scuff = squeeze(aaa(:,1))*1e9;
SP_scuff = zeros(Nport,Nport,length(ff_scuff));

for ii = 1:Nport
    rowshift = 2*Nport*(ii-1)+1;
    for jj = 1:Nport
        shift = rowshift + 2*(jj-1);
        SP_scuff(ii,jj,:) = aaa(:,shift+1) + 1j*aaa(:,shift+2);
    end
end
        
aaa = load('PLANAR_1coil_4port_zparam_scuff');
ZP_scuff = zeros(Nport,Nport,length(ff_scuff));

for ii = 1:Nport
    rowshift = 2*Nport*(ii-1)+1;
    for jj = 1:Nport
        shift = rowshift + 2*(jj-1);
        ZP_scuff(ii,jj,:) = aaa(:,shift+1) + 1j*aaa(:,shift+2);
    end
end

save('PLANAR_1coil_4port_scuff.mat', 'ff_scuff', 'SP_scuff', 'ZP_scuff');