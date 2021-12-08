%
%
%  This script plots the results (i.e. the surface displacement )for  
%  the programm StrikeSlipDef 
%
%% Parameters
clear all 
kilo = 1000;
k_p_mon = 30.4167;
cm = 100;

%% Variables 
ratio = 0.5;  %1 ratio of strike slip to tensile (i.e. 1 => pure strike slip; 0.5=> 50% Tensile, 50% strike slip)
%% location of piezo-stations
lat_pzs = 40.81305;lon_pzs = 27.77177;
lat_pzn = 40.81552; lon_pzn = 27.7776;

%% location of geodetic stations
lat_tekr = 40.958; long_tekr = 27.496;

%% read in source locations  
f = load('../res/fault_latlon.txt');
[n,m] = size(f);
f_lon = zeros([n,4]);
f_lat = zeros([n,4]);
fz = zeros([n,4]);
stk = zeros([n,4]);

bary_e = zeros([n,4]);
bary_w = zeros([n,4]);
for i = 1:n
    f_lon(i,:) = f(i,1:3:10);
    f_lat(i,:) = f(i,2:3:11);
    fz(i,:) = f(i,3:3:12)./kilo;
    stk(i,:) = f(i,13);
    bary_e(i,:) = f(i,14);
    bary_w(i,:) = f(i,15);

    f_lon_ave(i) = sum(f_lon(i,1:2))/2.; 
%     f_lat_ave(i) = sum(f_lat(i,1:2))/2.; 
end

%% read in calculatitons 
% read in displacment due to shear dislocation 
ss_file = '../res/obs_tekr_2000.dat';
[W2Es_loc,W2Es_dist,W2Etekr_ux_ss,W2Etekr_uy_ss,W2Etekr_uz_ss,W2E_dux_ss,W2E_duy_ss,W2E_duz_ss]= W2E_geodetic(ss_file,f_lon_ave,bary_e,bary_w);
[E2Ws_loc,E2Ws_dist,E2Wtekr_ux_ss,E2Wtekr_uy_ss,E2Wtekr_uz_ss,E2W_dux_ss,E2W_duy_ss,E2W_duz_ss]= E2W_geodetic(ss_file ,f_lon_ave,bary_e,bary_w);

% read in displacment due to tensile dislocation 
ten_file = '../res/obs_tekr_ten_2000.dat';
[W2Es_loc,W2Es_dist,W2Etekr_ux_ten,W2Etekr_uy_ten,W2Etekr_uz_ten,W2E_dux_ten,W2E_duy_ten,W2E_duz_ten]= W2E_geodetic(ten_file,f_lon_ave,bary_e,bary_w);
[E2Ws_loc,E2Ws_dist,E2Wtekr_ux_ten,E2Wtekr_uy_ten,E2Wtekr_uz_ten,E2W_dux_ten,E2W_duy_ten,E2W_duz_ten]= E2W_geodetic(ten_file,f_lon_ave,bary_e,bary_w);

% combine to produce final surface displacement at geodetic station 
W2E_ux = ratio.*W2Etekr_ux_ss+(1-ratio).*W2E_dux_ten;
W2E_uy = ratio.*W2Etekr_uy_ss+(1-ratio).*W2E_duy_ten;
W2E_uz = ratio.*W2Etekr_uz_ss+(1-ratio).*W2E_duz_ten;

E2W_ux = ratio.*E2Wtekr_ux_ss+(1-ratio).*E2W_dux_ten;
E2W_uy = ratio.*E2Wtekr_uy_ss+(1-ratio).*E2W_duy_ten;
E2W_uz = ratio.*E2Wtekr_uz_ss+(1-ratio).*E2W_duz_ten;


%% plot results 
 [n,m] =  size(W2Etekr_ux_ss);  % required for all plots 

figure(1)
set(gcf,'color','w');
hold on 
plot(W2Es_loc,W2E_ux.*cm,'.','Color','r','LineWidth',2)
plot(W2Es_loc,W2E_uy.*cm,'--','Color','b','LineWidth',2)
plot(W2Es_loc,W2E_uz.*cm,'-','Color','k','LineWidth',2)
% plot location of geodetic station and pressure sensors 
plot(long_tekr,0,'^k','MarkerSize',12,'LineWidth',2)
plot(lon_pzs,0,'.k','MarkerSize',14,'LineWidth',2)
plot(lon_pzn,0,'.k','MarkerSize',14,'LineWidth',2)
hold off 
legend('E-W','N-S disp.','Up')
title(['West to East, strike slip: ',num2str(ratio*100,'%3.1f'),'% tensile: ',num2str((1-ratio)*100,'%3.1f'),'%'])

xlabel('Source Longitude')
ylabel('Displacement (cm)')
grid on; box on 
pbaspect([2 1 1])

figure(2)
set(gcf,'color','w');
hold on 
plot(E2Ws_loc,E2W_ux.*cm,'.','Color','r','LineWidth',2)
plot(E2Ws_loc,E2W_uy.*cm,'--','Color','b','LineWidth',2)
plot(E2Ws_loc,E2W_uz.*cm,'-','Color','k','LineWidth',2)
% plot location of geodetic station and pressure sensors 
plot(long_tekr,0,'^k','MarkerSize',12,'LineWidth',2)
plot(lon_pzs,0,'.k','MarkerSize',14,'LineWidth',2)
plot(lon_pzn,0,'.k','MarkerSize',14,'LineWidth',2)
hold off 
legend('E-W','N-S disp.','Up')
title(['East to West, strike slip: ',num2str(ratio*100,'%3.1f'),'% tensile: ',num2str((1-ratio)*100,'%3.1f'),'%'])
xlabel('Source Longitude')
ylabel('Displacement (cm)')
grid on; box on 
pbaspect([2 1 1])


% 
figure(3)
set(gcf,'color','w');
hold on
plot(E2Ws_loc,E2W_uz./E2W_uy,'Color','k','LineWidth',2)
% plot location of geodetic station and pressure sensors 
plot(long_tekr,0,'^k','MarkerSize',12,'LineWidth',2)
plot(lon_pzs,0,'.k','MarkerSize',14,'LineWidth',2,'Color','b')
plot(lon_pzn,0,'.k','MarkerSize',14,'LineWidth',2,'Color','r')
hold off
title(['East to West, strike slip: ',num2str(ratio*100,'%3.1f'),'% tensile: ',num2str((1-ratio)*100,'%3.1f'),'%'])
xlabel('Source Longitude')
ylabel('Up/North Ratio')
set(gca, 'XDir','reverse')
axis([min(E2Ws_loc) max(E2Ws_loc) -3 3])
grid on; box on 



figure(4)
set(gcf,'color','w');
hold on 
vr = k_p_mon*0.53;  %   km/month
Prop_time = E2Ws_dist./kilo/vr;
plot(Prop_time,E2W_uy.*cm,'--','Color','k','LineWidth',2)
plot(Prop_time,E2W_uz.*cm,'-','Color','b','LineWidth',2)
grid on; box on
legend('North','Vertical')
xlabel('Time [months]')
ylabel('Displacement [cm]')
title(['East to West, strike slip: ',num2str(ratio*100,'%3.1f'),'% tensile: ',num2str((1-ratio)*100,'%3.1f'),'%'])




%% functions for reading files 
function[W2E_sources,dist_sources,W2E_ux,W2E_uy,W2E_uz,dux,duy,duz]= W2E_geodetic(fname,f_lon_ave,bary_x,bary_y)
%  Function calculating the surface displacement at point for a source moving from
%  West to East along the Marmara fault. 
%     fname:            file name for output from StrikeSlipDef programme
%     f_lon_ave:        corresponding longitude of sequence of sources
%                       along the Marmara fault 

% Outputs:
%     *_sources:        location of active source
%     *_sources:        distance in m from start of fault trace
%     *_ux,*_uy,*_uz:   permanent displacement at geodetic station due to accumlative slip along fault 
%     *du,*dy,*dz:      displacement at geodetic station at active cell   
%

a=fopen(fname);
rline = fgetl(a); 
n_sources = sscanf(rline, '%d');
rline = fgetl(a); 
n_obs = sscanf(rline, '%d');
rline = fgetl(a); 
dz = sscanf(rline, '%f');
rline = fgetl(a); 
data = sscanf(rline, '%f %f');
Lat = data(2);
Lon = data(1);
rline = fgetl(a); 
dux = sscanf(rline, '%f');
rline = fgetl(a); 
duy = sscanf(rline, '%f');
rline = fgetl(a); 
duz = sscanf(rline, '%f');
fclose(a);

E2Wtekr_dux = dux';
E2Wtekr_duy = duy';
E2Wtekr_duz = duz';


nn = length(tekr_duy);
W2E_sources = fliplr(f_lon_ave);
w2e_bary_x = fliplr(bary_x');
w2e_bary_y= fliplr(bary_y');

for i = 1:length(f_lon_ave)
   dist_sources(i) = sqrt((w2e_bary_x(i)-w2e_bary_x(1)).^2+(w2e_bary_y(i)-w2e_bary_y(1)).^2);
end


for i = 1:nn
   W2E_ux(i) = sum(E2Wtekr_dux(nn+1-i:nn));
   W2E_uy(i) = sum(E2Wtekr_duy(nn+1-i:nn));
   W2E_uz(i) = sum(E2Wtekr_duz(nn+1-i:nn));
end

end  


function[E2W_sources,dist_sources,E2W_ux,E2W_uy,E2W_uz,dux,duy,duz]= E2W_geodetic(fname,f_lon_ave,bary_x,bary_y)
%  Function calculating the surface displacement at point for a source moving from
%  East to West along the Marmara fault. 
%     fname:            file name for output from StrikeSlipDef programme
%     f_lon_ave:        corresponding longitude of sequence of sources
%                       along the Marmara fault 
% Outputs:
%     *_sources:        location of active source
%     *_sources:        distance in m from start of fault trace
%     *_ux,*_uy,*_uz:   permanent displacement at geodetic station due to accumlative slip along fault 
%     *du,*dy,*dz:      displacement at geodetic station at active cell   
%

a=fopen(fname);
rline = fgetl(a); 
n_sources = sscanf(rline, '%d');
rline = fgetl(a); 
n_obs = sscanf(rline, '%d');
rline = fgetl(a); 
dz = sscanf(rline, '%f');
rline = fgetl(a); 
data = sscanf(rline, '%f %f');
Lat = data(2);
Lon = data(1);
rline = fgetl(a); 
dux = sscanf(rline, '%f');
rline = fgetl(a); 
duy = sscanf(rline, '%f');
rline = fgetl(a); 
duz = sscanf(rline, '%f');
fclose(a);

E2Wtekr_dux = dux';
E2Wtekr_duy = duy';
E2Wtekr_duz = duz';

nn = length(E2Wtekr_duy);
for i = 1:length(f_lon_ave)
   dist(i) = sqrt((bary_x(i)-bary_x(1)).^2+(bary_y(i)-bary_y(1)).^2);
end

ii = 0;
clear dux duy duz
for i = 1:nn
   if (f_lon_ave(i) < 31.)
       ii = ii+1;
       E2W_sources(ii) = f_lon_ave(i);
       dist_sources(ii) = dist(i);
       dux(ii) = E2Wtekr_dux(i);
       duy(ii) = E2Wtekr_duy(i);
       duz(ii) = E2Wtekr_duz(i);
       if ii == 1
           E2W_ux(ii) = E2Wtekr_dux(i); % sum(E2Wtekr_dux(1:i));
           E2W_uy(ii) = E2Wtekr_duy(i); % sum(E2Wtekr_duy(1:i));
           E2W_uz(ii) = E2Wtekr_duz(i); % sum(E2Wtekr_duz(1:i));           
       else 
           E2W_ux(ii) = E2W_ux(ii-1)+E2Wtekr_dux(i); % sum(E2Wtekr_dux(1:i));
           E2W_uy(ii) = E2W_uy(ii-1)+E2Wtekr_duy(i); % sum(E2Wtekr_duy(1:i));
           E2W_uz(ii) = E2W_uz(ii-1)+E2Wtekr_duz(i); % sum(E2Wtekr_duz(1:i));
       end
   end
end

end  
