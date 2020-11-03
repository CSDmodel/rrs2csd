clear;

infile=dir('test.nc');

[num kk]=size(infile);
dc=360;dl=180;

aw412=0.00464;aw443=0.007098;aw469=0.010499;aw488=0.014579;
aw531=0.043935;aw547=0.053375;aw555=0.059649;aw645=0.326658;
aw667=0.433916;aw678=0.460971;bw412=0.003327;bw443=0.002438;
bw469=0.001909;bw488=0.001611;bw531=0.001123;bw547=0.000989;
bw555=0.00093;bw645=0.000490;bw667=0.000425;bw678=0.000397;

po412=555.0/412.0;po443=555.0/443.0;po469=555.0/469.0;po488=555.0/488.0;
po531=555.0/531.0;po547=555.0/547.0;po555=555.0/555.0;po645=555.0/645.0;
po667=555.0/667.0;po678=555.0/678.0;

lf1=[-0.209,0.426,-0.567,0.312,-0.185,0.523,-0.092,-0.418,0.064]; 
lf2=[-0.285,0.215,-0.457,0.404,-0.128,0.300,-0.059,-0.275,0.100];

for k=1:num
    input=char(infile(k).name)
    
    rrs412=double(ncread(input,'Rrs_412'));
    rrs443=double(ncread(input,'Rrs_443'));
    rrs469=double(ncread(input,'Rrs_469'));
    rrs488=double(ncread(input,'Rrs_488'));
    rrs531=double(ncread(input,'Rrs_531'));
    rrs547=double(ncread(input,'Rrs_547'));
    rrs555=double(ncread(input,'Rrs_555'));
    rrs667=double(ncread(input,'Rrs_667'));
    lat=double(ncread(input,'lat'));
    lon=double(ncread(input,'lon'));
    
    br=log10(max(rrs443,rrs488)./rrs547);
    chl=10.^(0.2424-2.7423*br+1.8017*(br.^2)+0.0015*(br.^3)-1.2280*(br.^4)); 
    
    ulim667=20.*(rrs555).^1.5;llim667=0.9.*(rrs555).^(1.7);
  
    if (rrs667>=ulim667)
        srrs667=ulim667./(0.52+1.7*ulim667);
    elseif (rrs667<=llim667)
        srrs667=llim667./(0.52+1.7*llim667);
    else
        srrs667=rrs667./(0.52+1.7*rrs667);
    end
    
    srrs412=rrs412./(0.52+1.7*rrs412);
    srrs443=rrs443./(0.52+1.7*rrs443);
    srrs469=rrs469./(0.52+1.7*rrs469);
    srrs488=rrs488./(0.52+1.7*rrs488);
    srrs531=rrs531./(0.52+1.7*rrs531);
    srrs547=rrs547./(0.52+1.7*rrs547);
    srrs555=rrs555./(0.52+1.7*rrs555);
    
    g0=0.089;g1=0.125;
    
    u412=(-g0+sqrt(g0^2+4*g1.*srrs412))/(2*g1);
    u443=(-g0+sqrt(g0^2+4*g1.*srrs443))/(2*g1);
    u469=(-g0+sqrt(g0^2+4*g1.*srrs469))/(2*g1);
    u488=(-g0+sqrt(g0^2+4*g1.*srrs488))/(2*g1);
    u531=(-g0+sqrt(g0^2+4*g1.*srrs531))/(2*g1);
    u547=(-g0+sqrt(g0^2+4*g1.*srrs547))/(2*g1);
    u555=(-g0+sqrt(g0^2+4*g1.*srrs555))/(2*g1);
    u667=(-g0+sqrt(g0^2+4*g1.*srrs667))/(2*g1);
          
    X=log10((srrs443+srrs488)./(srrs555+(5*(srrs667.^2)./srrs488)));
    Y=2.0*(1.0-1.2*exp(-0.9*srrs443./srrs555));
    Z=0.74+0.2*(0.8+srrs443./srrs555).^(-1);
    S=0.015+0.002*(0.6+(srrs443./srrs555)).^(-1);
    E=exp(S*(443.0-412.0));
    a555=aw555+10.^(-1.146-1.366*X-0.469*X.^2);
    bp555=(u555.*a555)./(1.0-u555)-bw555;
    
    bp412=bp555.*(po412).^Y;bp443=bp555.*(po443).^Y;bp469=bp555.*(po469).^Y;
    bp488=bp555.*(po488).^Y;bp531=bp555.*(po531).^Y;bp547=bp555.*(po547).^Y;
    bp667=bp555.*(po667).^Y;
    b412=bw412+bp412;b443=bw443+bp443;b469=bw469+bp469;b488=bw488+bp488;
    b531=bw531+bp531;b547=bw547+bp547;b555=bw555+bp555;b667=bw667+bp667;
    
    a412=(1.0-u412).*b412./u412;a443=(1.0-u443).*b443./u443;
    a469=(1.0-u469).*b469./u469;a488=(1.0-u488).*b488./u488;
    a531=(1.0-u531).*b531./u531;a547=(1.0-u547).*b547./u547;
    a555=(1.0-u555).*b555./u555;a667=(1.0-u667).*b667./u667;
    
    adg443=((a412-Z.*a443)-(aw412-Z.*aw443))./(E-Z);
    
    aph412=real(a412-aw412-adg443.*exp(-S*(412-443)));
    aph443=real(a443-aw443-adg443);
    aph469=real(a469-aw469-adg443.*exp(-S*(469-443)));
    aph488=real(a488-aw488-adg443.*exp(-S*(488-443)));
    aph531=real(a531-aw531-adg443.*exp(-S*(531-443)));
    aph547=real(a547-aw547-adg443.*exp(-S*(547-443)));
    aph555=real(a555-aw555-adg443.*exp(-S*(555-443)));
    aph667=real(a667-aw667-adg443.*exp(-S*(667-443)));
    
    data(:,:,1)=aph412;data(:,:,2)=aph443;data(:,:,3)=aph469;data(:,:,4)=aph488;
    data(:,:,5)=aph531;data(:,:,6)=aph547;data(:,:,7)=aph555;
    
    aph_std = std(data,0,3);
    aph_mean = mean(data,3);

    data2 = NaN(length(lon),length(lat),8);
    csd = NaN(length(lon),length(lat));
    for iLon = 1 : length(lon)
        for iLat = 1 : length(lat)
            if  (aph412(iLon,iLat) > aph469(iLon,iLat))
                data2(iLon,iLat,1)=((data(iLon,iLat,1)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(2);
                data2(iLon,iLat,2)=((data(iLon,iLat,2)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(3);
                data2(iLon,iLat,3)=((data(iLon,iLat,3)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(4);
                data2(iLon,iLat,4)=((data(iLon,iLat,4)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(5);
                data2(iLon,iLat,5)=((data(iLon,iLat,5)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(6);
                data2(iLon,iLat,6)=((data(iLon,iLat,6)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(7);
                data2(iLon,iLat,7)=((data(iLon,iLat,7)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf1(8);
                data2(iLon,iLat,8)=chl(iLon,iLat).*lf1(9);
                csd(iLon,iLat)=(lf1(1)+exp(sum(data2(iLon,iLat,:)))).^(-1);
            elseif (aph412(iLon,iLat) < aph469(iLon,iLat))
                data2(iLon,iLat,1)=((data(iLon,iLat,1)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(2);
                data2(iLon,iLat,2)=((data(iLon,iLat,2)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(3);
                data2(iLon,iLat,3)=((data(iLon,iLat,3)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(4);
                data2(iLon,iLat,4)=((data(iLon,iLat,4)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(5);
                data2(iLon,iLat,5)=((data(iLon,iLat,5)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(6);
                data2(iLon,iLat,6)=((data(iLon,iLat,6)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(7);
                data2(iLon,iLat,7)=((data(iLon,iLat,7)-aph_mean(iLon,iLat))./aph_std(iLon,iLat)).*lf2(8);
                data2(iLon,iLat,8)=chl(iLon,iLat).*lf2(9);
                csd(iLon,iLat)=(lf2(1)+exp(sum(data2(iLon,iLat,:)))).^(-1);
            end
        end
    end

    ind=find(csd>=3.0);csd(ind)=3;
    ind=find(csd<=0.0);csd(ind)=0;
end
