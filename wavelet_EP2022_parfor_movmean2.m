
%uses croscabernet2021_lod, with smoothing and local maximum for FFT, and
%flux LOD
%Dec 10, 2021 added moving mean
% temporary change for flight 1: altagl
%Jan 24, 2022: added new random flux calc
% change Feb 4, 2022: changed dt calculation from mean to mode (most common
% time interval instead of avg) & added flag for if covariance peak is just
% noise
% Feb23,2022: changed LOD calculation and discarded endpoints of moving
% mean
% Mar 4, 2022: changed the code so that the flux calculation is done using
% ppb - conversion to mass flux afterwards

% 1)----load and prepare data--------------------------%

addpath(genpath('C:\Users\evayp\OneDrive\Documents\Matlab'))
clear

segx=importdata('segmentindex.mat');

numsegm =length(segx)/2; % number of segments


parfor segm=1:numsegm
segx=importdata('segmentindex.mat');
numsegm =length(segx)/2;
data_clean=importdata('data_clean.mat');
cir=load('cirpas.mat');
cirpas = cir.cirpas;                   
cirpas_i = cir.cirpas_i;               
gpstimedoy = cir.gpstimedoy;           
gpstimedoy_LT = cir.gpstimedoy_LT;
pblinterpolated=importdata('pblinterpolated.mat','pblinterpolated');
if isstruct(pblinterpolated)
    if isfield(pblinterpolated,'idxoutpbl')
    idxoutpbl = pblinterpolated.idxoutpbl;
    else
    idxoutpbl = zeros(length(pblinterpolated.pblinterpolated),1)
    end
    pblinterpolated =pblinterpolated.pblinterpolated;
elseif ~isstruct(pblinterpolated)
    idxoutpbl = zeros(length(pblinterpolated),1)
end


w_star=ones(length(data_clean.UTCTime),1).*NaN; %empty cells for w* = convective velocity scale = Deardorff velocity
O=struct;
O.ciri(:,1)=data_clean.Lat; %aircraft data, Lat, interpolated on ptr time
O.ciri(:,2)=data_clean.Long;%aircraft data, Lon, interpolated on ptr time
O.RadT_C=data_clean.IRTNad_C_; %aircraft data, temperature, interpolated on ptr time
%O.pblinterpolatedv=interp1(cirdoy,pblinterpolated,ptrdoy); %aircraft data, PBL, interpolated on ptr time
O.NovAtelAltmv=data_clean.NovAtelAlt_m_; %aircraft data, altitude, interpolated on ptr time
O.TASmsv=data_clean.TAS_m_s_; %aircraft data, true airspeed, interpolated on ptr time
O.RadAltm = data_clean.RadAlt__m_;
tstimeLTdoy=dec_doy(data_clean.UTCTime)-7/24; % bring on local time

%% 2)---determine altitude above ground ---------%
%scatter(O.RadAltm,O.NovAtelAltmv);
%plot(tstimeLTdoy,O.RadAltm,tstimeLTdoy,O.NovAtelAltmv);
%radar altitude doesn't work above 700 m (?), so this is padded with data from
%sensor that measures height a.s.l.
O.RadAltmValid=filtfilt(ones(200,1)/200,1,O.RadAltm); %filtfilt function smoothes data without shifting vs x
O.RadAltmValid(O.RadAltmValid>700)=NaN;
 rr=isnan(O.RadAltmValid);
 doynotnan=tstimeLTdoy(rr==0);
 RadAltmValidnotnan=O.RadAltmValid(rr==0);
 NovAtelAltmnotnan=O.NovAtelAltmv(rr==0);
 Groundlevel=NovAtelAltmnotnan-RadAltmValidnotnan;
 Groundlevel=interp1(doynotnan,Groundlevel,tstimeLTdoy);
O.Altagl=O.NovAtelAltmv-Groundlevel;
% alg = importdata("altagl.mat");
% O.Altagl = alg.Altagl;

VertWindms = cirpas.Vert_Wind_m_s_;
VertWindmsidx =(VertWindms == -9999  | VertWindms > 50);
VertWindms(VertWindmsidx) = NaN;
 

 %% Calculate mg/m³ from ppb, and sampling intervals

amus = readtable('amus.txt'); % get m/z values from txt exported from Tofware
Peaknames = readtable('PeakTable.csv');
amus_noProt = amus.amus-1.008;
names = data_clean(:,1:end-32).Properties.VariableNames;
O.presi=data_clean.Ps_mb_;%interpolated pressure from aircraft data to PTR time
O.tempi=data_clean.Tamb_C_;%interpolated T from aircraft data to PTR time

VOCs_ppb = table2array(data_clean(:,1:end-32));
%data_clean_mgm3 = (table2array(data_clean(:,1:end-32))); %.*(amus_noProt'))./(1013*22.4*(273+O.tempi(1:size(O.tempi)))./(273*O.presi(1:size(O.presi))))*3600./1000;
% do the flux calculation with ppb*m/s, don't want to introduce any bias
% from the temperature fluctuations
%data_clean_ppb(end,:)=0;

O.dt=mode(diff(gpstimedoy_LT*3600*24)); %most frequent sampling interval cirpas
O.dt_voc=mode(diff(tstimeLTdoy*3600*24)); %most frequent sampling interval VOC

O.vwi=data_clean.Vert_Wind_m_s_; %interpolated vertical ws to ptr data

amus=importdata("amus_do.mat");
% LOD=load('LOD.mat');
% amus.do = zeros(length(table2array(amus)),1);

% FlagsPerCol = sum(LODflagmatrix,1)./size(LODflagmatrix,1);
% LODflagmass = FlagsPerCol >= 0.999;
% remainMass = length(amus_noProt) -sum(LODflagmass)
%amus.do(LODflagmass(1:end-32) ==0) = 1;
%save amus_do_LOD.mat amus
dofilter=(amus.do ==1 & amus.amus >30);
%dofilter=(amus.amus >30); % & amus.amus< 69.07);


namesVOC_do= names(dofilter);
amus_do = amus_noProt(dofilter);

VOCppb_do = VOCs_ppb(:,dofilter);
nVOC = size(VOCppb_do,2);

O.VOCmgm3i=interp1(tstimeLTdoy,VOCppb_do,gpstimedoy_LT); %interpolate vocs to aircraft time axis
%% loop through all masses

for v=1:nVOC

O.VOCppb = VOCppb_do(:,v);
O.amu = amus_do(v);



 %% Remove spikes in roll 
%1) for data in PTR time resolution

rollmax = 8; %pitchmin = 4
nanwindow=3;
O.RollPTR= data_clean.Roll_deg_;
O.RollIdxPTR = (abs(data_clean.Roll_deg_) >= rollmax);
O.PTRds = find(abs(data_clean.Roll_deg_) >= rollmax) ;

for i=-nanwindow:nanwindow
    if (O.PTRds+i)<=length(O.VOCppb) & (O.PTRds+i) >0
     O.VOCppb(O.PTRds+i) = NaN;
    O.RollPTR(O.PTRds+i)=NaN;
    end
end 

%make a stitched version of VOC and its time axis for filted plot
ptrdoynotnan = tstimeLTdoy;
O.VOCppbnotnan = O.VOCppb;
bad = isnan(O.VOCppb);
O.VOCppbnotnan(bad)=[]; ptrdoynotnan(bad)=[];
 
% 2) for data in cirpas time resolution
O.cirds = find(abs(cirpas.Roll_deg_) >= rollmax);

for i=-(nanwindow):(nanwindow)
    if (O.cirds+i) <= length(gpstimedoy_LT) & (O.cirds+i) >0
    VertWindms(O.cirds+i) = NaN;
    %TambCOrig(O.cirds+i) = NaN;
    end
end



%% make subfolders to store figures
D = cd;
F = 'Segments';
if ~exist(fullfile(D,F),'dir')
    mkdir(fullfile(D,F))
end

    Z = fullfile(D,F,sprintf('Segment%d',segm));
    if ~exist(Z,'dir')
        mkdir(Z)
    end
    
     Y = fullfile(Z,sprintf(namesVOC_do{v}));
    if ~exist(Y,'dir')
        mkdir(Y)
    end
%% -------select flight segment---------------------

%for PBL height determination set pblflagt to 1, select start and end of
%each sounding, place marker on the point in the sounding which corresponds
%to PBL height. The two such markers will be used for interpolation

%[xx,yy]=ginput(16);
%    xx(1) = tstimeLTdoy(indptr(1));
%  xx(2)=tstimeLTdoy(indptr(end));
%save segmentindex.mat xx

%segx=importdata('segmentindex.mat');
xx=[];

 xx(1) = segx(segm*2-1);
 xx(2)=segx(segm*2);

% klick markers to determine beginning and end of flight stretch to
% analyze. Should not contain turns or significant changes in roll/altitude
O.indcir=find(gpstimedoy_LT(:)>xx(1) & gpstimedoy_LT(:)<xx(2));%determines index for start and end of chosen segment
O.indptr=find(tstimeLTdoy(:)>xx(1) & tstimeLTdoy(:)<xx(2)); %determines index for start and end of chosen segment
O.segmentlengthkm =(xx(2)-xx(1))*24*60*60*mean(cirpas.TAS_m_s_(O.indcir))/1000;
O.segmenttimeminutes =(xx(2)-xx(1))*24*60 ;


% save segment data in O for future reference:
O.segm(1,1)=O.segmentlengthkm; O.segm(1,2)=O.segmenttimeminutes; O.segm(1,3)=gpstimedoy_LT(O.indcir(1)); O.segm(1,4)=gpstimedoy_LT(O.indcir(end));O.segm(1,5)=data_clean.Lat(O.indptr(1)); O.segm(1,6)=data_clean.Long(O.indptr(1)); O.segm(1,7)=data_clean.Lat(O.indptr(end)); O.segm(1,8)=data_clean.Long(O.indptr(end));


%% ----export data from selected segment into following variables for wavelet/flux analysis-----#
% excluding the NaNs from roll despiking

O.VOC=[];
O.vw=[];
O.idxSegment = ~isnan(O.VOCppb(O.indptr));
O.VOC(:,1)=tstimeLTdoy(O.indptr(O.idxSegment)); %1st column in O.VOC is ptr time
O.VOC(:,2)=O.VOCppb(O.indptr(O.idxSegment)); %2nd column is VOC in ppb
O.idxSegmentvw = ~isnan(VertWindms(O.indcir));
O.vw(:,1)=gpstimedoy_LT(O.indcir(O.idxSegmentvw)); %1st column of vertical windsp is aircraft time
O.vw(:,2)=VertWindms(O.indcir(O.idxSegmentvw)); %2nd col is vertical windsp
O.vw(:,3)=cirpas.Theta_K_(O.indcir(O.idxSegmentvw)); %3rd col is potential temperature


%% ---wavelet transform-----------------------------
O.pw=O.dt_voc; %dt is (disjunct) sampling time, i.e. data interval
O.s_0_voc=O.dt_voc*2; % 2 x ptrms datainterval= scaling factor
O.s_0=0.2; %2 x 0.1 s for 10 Hz data
O.dj=0.25; %octaves
O.delay=3.3;
O.flag=1; %plot wavelet figure in function? 
O.Ts=4; %integral timescale: s. integral timescale as Gamma in Karl 2013, appendix B, is only used to calc. disjunct error
O.TASkms=mean(O.TASmsv(O.indptr)/1000); %aircraft speed: km/s
O.pbl=mean(pblinterpolated(O.indptr)/1000); %1.5; %km 
O.zi=mean(O.Altagl(O.indptr)/1000); %km

%wavecroscabernet calculates fft and wavelet 
spec = wavecroscabernet_lod2021(O.VOC, O.vw, O.delay, O.s_0, O.s_0_voc, O.dj, O.zi, O.pbl, O.TASkms, O.Ts,O.flag);

exportgraphics(gcf,(fullfile(Y,'covpeak.png')),'Resolution',300)
%saveas(gcf,(fullfile(Y,'covpeak.png')))  
close all
% % convective velocity scale = deardorff velocity is calculated in
% % wavecreoscabernet and now interpolated on ptr spectra time
% wstarptr=interp1(spec.timew,spec.wstar,spec.timx)';
% w_star(O.indptr(1):O.indptr(1)+length(wstarptr)-1)=wstarptr; 
% %saved in w_star but only for the segment that was chosen

    
%% ----plot comparison between FFT and wavelet transform and model--------


% ---1) model: see Handbook of Micrometeorology, p.45---
% a0: normalization parameter
% m: intertial subrange slope parameter
% mi: broadness parameter. smaller mi = narrower cospectra, approaching
% flat-terrain stable-atmosphere; the opposite for larger mi
% fx: variability in the frequency, fx, at which the frequencyweighted 
% cospectra, fCo(f), reaches a maximum value

a0=0.4;%3
m=1/8;%1/2 %1/3
mi=4/3;%7/3 %5/3
fx=0.009;%0.05  %0.025
modelcos=a0.*((spec.fqfft./fx)./(1+m.*spec.fqfft./fx).^(2*mi)).^((1/(2*mi))*((m+1)/m));
nyquist = 1/(O.dt_voc*2); %Nyquist freq of VOC measurement (half the sampling freq)


% 2) plot normalized (to maximum) co-spectra of VOC and heat flux. Wavelet
% spectra are multiplied with frequency to remove frequency dependence
% included in y (?) and get to 2d

ax=figure; semilogx(spec.fqwavevoc,spec.fqwavevoc'.*spec.wavewVOC./size(spec.wavewVOC,1)./...
     max(spec.fqwavevoc'.*spec.wavewVOC./size(spec.wavewVOC,1)),'linewidth',3);
     hold on; plot(ones(22,1)*nyquist,[-1:0.1:1.1]),'k'; semilogx(spec.fqfft,modelcos,'--','linewidth',1);
 xlabel('Hz'); ylabel('f x CoS'); legend('<w'' ''VOC'' ''_w_a_v_e>','Nyquist frequency','model','Location','southeast')
C = colororder(ax);
C(3,:) = C(1,:); C(4,:)=C(2,:);C(5,:)=[0 0 0];
colororder(ax,C);
% ax=figure; semilogx(spec.fqfft,spec.wT./max(spec.wT),'--',spec.fqfft, spec.wvoc./max(spec.wvoc),'--',...
%      spec.fqwave,spec.fqwave'.*spec.wavewT./max(spec.fqwave'.*spec.wavewT),spec.fqwavevoc,spec.fqwavevoc'.*spec.wavewVOC./size(spec.wavewVOC,1)./...
%      max(spec.fqwavevoc'.*spec.wavewVOC./size(spec.wavewVOC,1)),'linewidth',3);
%      hold on; plot(ones(22,1)*nyquist,[-1:0.1:1.1]),'k'; semilogx(spec.fqfft,modelcos,'--','linewidth',1);
%  xlabel('Hz'); ylabel('f x CoS'); legend('<w'' ''T'' ''_f_f_t>', '<w'' ''VOC'' ''_f_f_t>', '<w'' ''T'' ''_w_a_v_e>','<w'' ''VOC'' ''_w_a_v_e>','Nyquist frequency','model','Location','southeast')
% C = colororder(ax);
% C(3,:) = C(1,:); C(4,:)=C(2,:);C(5,:)=[0 0 0];
% colororder(ax,C);

% plot second x axis in km (reverse), linked to scale of frequency axis
rr1=floor(log10((min(spec.fqfft(spec.fqfft>0)))));
ll1=ceil(log10((max(spec.fqfft))));
LinkTopAxisData (10.^(rr1:ll1),O.TASkms./10.^(rr1:ll1),'km');
grid on

saveas(gcf,(fullfile(Y,['cos.png'])))  


% 3) plot cumulative co-spectra for VOC and heat fluxes
% yy1-4 are normalized to maximum and to the number of datapoints because
% FFT results have less datapoints than wavelet

%fqwave: wavelet frequency
%wavewT: (variancew*varianceT)*sum(wc_wT')/nw;
yy1=spec.fqwave'.*spec.wavewT./max(spec.fqwave'.*spec.wavewT);
yy2=spec.fqwavevoc'.*spec.wavewVOC./size(spec.wavewVOC,1)./...
     max(spec.fqwavevoc'.*spec.wavewVOC./size(spec.wavewVOC,1));
yy3=cumsum(spec.wT(~isnan(spec.wT))./max(spec.wT))/size(    spec.fqfft(~isnan(spec.wT)),1);
yy4=cumsum(spec.wvoc(~isnan(spec.wvoc))./max(spec.wvoc))/size(spec.fqfft(~isnan(spec.wvoc)),1);
yy5=cumsum(modelcos./max(spec.wvoc))/size(spec.fqfft(~isnan(spec.wvoc)),1);
yy1=cumsum(flipud(yy1))/size(yy1,1);
yy2=cumsum(flipud(yy2))/size(yy2,1);

 figure; semilogx(fliplr(spec.fqwavevoc),yy2/yy2(size(yy2,1)),'linewidth',3);
     hold on; plot(ones(22,1)*nyquist,[-1:0.1:1.1]);semilogx(spec.fqfft, yy5/yy5(size(yy5,1)),'k--')
 xlabel('Hz'); ylabel('\int f x CoS'); legend('<w'' ''VOC'' ''_w_a_v_e>','Nyquist frequency','model',"location","southwest")
 grid on

%  figure; semilogx(fliplr(spec.fqwave),yy1/yy1(size(yy1,1)),fliplr(spec.fqwavevoc),yy2/yy2(size(yy2,1)),'linewidth',3);
%      hold on; plot(ones(22,1)*nyquist,[-1:0.1:1.1]);semilogx(spec.fqfft, yy5/yy5(size(yy5,1)),'k--')
%  xlabel('Hz'); ylabel('\int f x CoS'); legend('<w'' ''T'' ''_w_a_v_e>','<w'' ''VOC'' ''_w_a_v_e>','Nyquist frequency','model',"location","southwest")
%  grid on

 saveas(gcf,(fullfile(Y,['cumcos.png'])))  
%% plot flux (i.e. integral over all frequency scales)
coimax = 0.8;

O.fftwavratio=spec.fftfluxvoc./mean(spec.wavefluxvoc);
 

%Now attention: spec has less datapoints than O.Spr due to NANs - so need
%to match timepoints of other parameters on the correct times
spec.timxd=spec.timx./(3600*24); %time difference
O.Spr = [];
 O.Spr(:,1)=spec.timxd+O.segm(1,3); %(O.segm(1,3) is the start time of the segment)
 O.Spr(:,2)=interp1(tstimeLTdoy,O.VOCppb,O.Spr(:,1),'nearest'); 
 O.Spr(:,3:4)=interp1(tstimeLTdoy,O.ciri(:,1:2),O.Spr(:,1),'nearest'); 
 O.Spr(:,5)=spec.fftfluxvoc;
 O.Spr(:,6)=spec.wavefluxvoc; 
 O.Spr(:,8)=interp1(tstimeLTdoy,O.NovAtelAltmv,O.Spr(:,1),'nearest') ;
 O.Spr(:,9)=interp1(tstimeLTdoy,O.Altagl,O.Spr(:,1)) ;
 O.Spr(:,10)=interp1(tstimeLTdoy,pblinterpolated,O.Spr(:,1),'nearest') ;

 O.Spr(:,15)=spec.fluxSNRwave ; %signal to noise ratio
 O.Spr(:,16)=spec.wavelet_noise_std ; %random flux std
 LODfluxindex =abs(spec.wavefluxvoc) <= spec.wavelet_noise_std;
 O.Spr(:,17)=LODfluxindex;
 SNRindex = spec.fluxSNRwave <=1;

figure; plot(spec.timx,spec.wavefluxvoc); hold on
plot(spec.timx,ones(size(spec.timx,2),1)*spec.fftfluxvoc)
plot(spec.timx,ones(size(spec.timx,2),1)*spec.wavelet_noise_std*1,'k',spec.timx,ones(size(spec.timx,2),1)*spec.wavelet_noise_std*(-1),'k')
plot(spec.timx,ones(size(spec.timx,2),1)*spec.random_flux*3,'g',spec.timx,ones(size(spec.timx,2),1)*spec.random_flux*(-3),'g')
ylabel('ppb*m/s')
legend('Wavelet flux','FFT flux','flux noise (Spirig 2005)','','3x FFT random flux (Billesbach 2011)','Location','northoutside')
%saveas(gcf,(fullfile(Y,['flux.png'])))


ziz=(O.Spr(:,9)./O.Spr(:,10));

 %equation to extrapolate flux to the ground (for MeOH no gradient assumed):
 isodiverg=1; %(ziz-1.1508)/(-1.1424); 
 spec.surfaceflux=spec.wavefluxvoc./isodiverg;
 spec.surfaceflux(abs(spec.qcoi)>coimax) = NaN;
 O.Spr(:,7)=spec.surfaceflux;
 figure; plot(O.Spr(:,1),O.Spr(:,6),O.Spr(:,1),O.Spr(:,7),'r',O.Spr(~LODfluxindex,1),O.Spr(~LODfluxindex,7)); legend('flux at z','surface flux','flux above LOD and COI','Location','best');
 
 O.Spr(:,11)=ziz;
 
 O.Spr(:,12) =spec.noiseflag;
 O.Spr(:,13)=interp1(tstimeLTdoy,O.RadT_C,O.Spr(:,1),'nearest') ; 
 O.Spr(:,14)=interp1(tstimeLTdoy,O.TASmsv,O.Spr(:,1),'nearest') ;

 % convective velocity scale = deardorff velocity is calculated in
% wavecreoscabernet and now interpolated on ptr spectra time
wstarptr=interp1(spec.timew,spec.wstar,spec.timx)';
O.Spr(:,18) = wstarptr;
O.Spr(:,19) = spec.se;
O.Spr(:,20) = spec.re;
O.Spr(:,21) = interp1(tstimeLTdoy,O.presi,O.Spr(:,1),'nearest') ; 
O.Spr(:,22)=interp1(gpstimedoy_LT,cirpas.WindSpeed_m_s_,O.Spr(:,1),'nearest') ;
O.Spr(:,23)=interp1(gpstimedoy_LT(O.indcir),cirpas.WindDir_Deg_(O.indcir),O.Spr(:,1),'nearest') ;%using "nearest" is ok for degree
O.Spr(:,24)=spec.lagtime ;
O.Spr(:,25)=interp1(spec.timew,spec.wavefluxwT,spec.timx)';
O.Spr(:,26)=interp1(tstimeLTdoy,O.tempi,O.Spr(:,1),'nearest') ; 
  
O.Spr(ziz>1,:) = []; % PBL filter especially important for LA
% interpolate manual pbl filter
if isvarname("idxoutpbl") 
    pblIndex = interp1(tstimeLTdoy,double(idxoutpbl),O.Spr(:,1));
    O.Spr(pblIndex==1,:) = [];
end


parsave(fullfile(Y,"O.mat"),O.Spr);



%% 2 km running avgs every 200 m
%should change this to filtfilt and downsample?
%spec.surfacefluxf = filtfilt(ones(5,1)/5,1,spec.surfaceflux)
O.t = 2000/mean(data_clean.TAS_m_s_(O.indptr))/3600/24;
O.tsub = 200/mean(data_clean.TAS_m_s_(O.indptr))/3600/24;
O.numavg = round(O.t/median(diff(O.Spr(:,1))));
O.numavgsub = round(O.tsub/median(diff(O.Spr(:,1))));
Out_doy = downsample(movmean(O.Spr(:,1),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_flux = downsample(movmean(O.Spr(:,7),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_ppb = downsample(movmean(O.Spr(:,2),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_TAS = downsample(movmean(O.Spr(:,14),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_PBLinterp = downsample(movmean(O.Spr(:,10),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_ziz = downsample(movmean(O.Spr(:,11),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_RadTC = downsample(movmean(O.Spr(:,13),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_wstar = downsample(movmean(O.Spr(:,18),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_windsp = downsample(movmean(O.Spr(:,22),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_lagtime= ones(length(Out_doy),1).*spec.lagtime;
Out_waveletnoise = downsample(movmean(O.Spr(:,16),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_wTflux = downsample(movmean(O.Spr(:,25),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_SE= ones(length(Out_doy),1).*spec.se;
Out_RE= ones(length(Out_doy),1).*spec.re;
Out_Tamb =  downsample(movmean(O.Spr(:,26),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
Out_press =  downsample(movmean(O.Spr(:,21),O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));

winddir= O.Spr(:,23);
    [x,y] = pol2cart(deg2rad(winddir),ones(size(winddir)));
    x=downsample(movmean(x,O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
    y=downsample(movmean(y,O.numavg,'omitnan'),O.numavgsub,floor(O.numavgsub/2));
    [windir_avged,~]=cart2pol(x,y);
    windir_avged = rad2deg(windir_avged);

Out_winddir = windir_avged;
% Out_fluxLODflag = abs(Out_flux) <= spec.wavelet_noise_std;
% Out_fluxLOD = Out_flux;
% Out_fluxLOD(Out_fluxLODflag) = NaN;

Out_Lat = interp1(tstimeLTdoy,data_clean.Lat, Out_doy,'nearest');
Out_Lon = interp1(tstimeLTdoy,data_clean.Long, Out_doy,'nearest');
Out_altagl = interp1(tstimeLTdoy,O.Altagl,Out_doy,'nearest');
SegNo = segm.*ones(length(Out_doy),1);

if v==1
VOCfluxes = zeros(length(Out_doy),nVOC);
VOCrandfluxes = zeros(length(Out_doy),nVOC);
VOCrandfluxesfft = zeros(length(Out_doy),nVOC);
VOCrandfluxesmov = zeros(length(Out_doy),nVOC);
Cov_noiseflag = zeros(length(Out_doy),nVOC);
VOC_lagtimes = zeros(length(Out_doy),nVOC);
VOC_ppb = zeros(length(Out_doy),nVOC);
ParSeg =table(Out_doy,Out_Lat,Out_Lon,Out_altagl,Out_TAS, Out_ziz, Out_PBLinterp,Out_RadTC,Out_wstar,Out_windsp,Out_winddir,Out_wTflux,Out_SE,Out_RE,SegNo);
end

if length(Out_flux) ~= size(VOCfluxes,1)
    %ParSeg =table(Out_doy,Out_Lat,Out_Lon,Out_altagl,Out_TAS, Out_ziz, Out_PBLinterp,Out_RadTC,Out_wstar,Out_windsp,Out_winddir,Out_wTflux,SegNo);
    Out_flux = interp1(Out_doy, Out_flux,table2array(ParSeg(:,1)));
    Out_press = interp1(Out_doy, Out_press,table2array(ParSeg(:,1)));
    Out_Tamb = interp1(Out_doy, Out_Tamb,table2array(ParSeg(:,1)));
end

    if length(Out_ppb) ~= size(VOCfluxes,1)
        Out_ppb = interp1(Out_doy, Out_ppb,table2array(ParSeg(:,1)));
    end

VOCfluxes(:,v) = Out_flux .*O.amu.*3600.*1000.*Out_press*100.*1e-9./(8.3144.*(237.15+Out_Tamb));
VOCfluxes_t=array2table(VOCfluxes);
VOCfluxes_t.Properties.VariableNames = namesVOC_do;
VOCrandfluxes(:,v) = spec.waveflux_sigma.*O.amu.*3600.*1000.*Out_press*100.*1e-9./(8.3144.*(237.15+Out_Tamb));
VOCrandfluxes_t=array2table(VOCrandfluxes);
VOCrandfluxes_t.Properties.VariableNames = append('sigma_',namesVOC_do);
VOCrandfluxesmov(:,v) = spec.wavelet_noise_movstd.*O.amu.*3600.*1000.*Out_press*100.*1e-9./(8.3144.*(237.15+Out_Tamb));
VOCrandfluxesmov_t=array2table(VOCrandfluxesmov);
VOCrandfluxesmov_t.Properties.VariableNames = append('sigmawav_',namesVOC_do);
Cov_noiseflag(:,v) = spec.noiseflag;
Cov_noiseflag_t=array2table(Cov_noiseflag);
Cov_noiseflag_t.Properties.VariableNames = append('Cov_noise_flag',namesVOC_do);
VOC_lagtimes(:,v) = spec.lagtime;
VOC_lagtimes_t = array2table(VOC_lagtimes);
VOC_lagtimes_t.Properties.VariableNames = append('Lagt_',namesVOC_do);


VOC_ppb(:,v) = Out_ppb;
VOC_ppb_t = array2table(VOC_ppb);
VOC_ppb_t.Properties.VariableNames = append('ppb_',namesVOC_do);

VOCSeg = horzcat(ParSeg,VOCfluxes_t,VOCrandfluxes_t,VOCrandfluxesmov_t,Cov_noiseflag_t,VOC_lagtimes_t,VOC_ppb_t);

close all
O.Spr =[];
end

parsave(sprintf('Seg%d.mat',segm),VOCSeg);
 
clc
   
end 

 %% After all segments have been processed
for k = 1:numsegm 
  Segs{k} = load(sprintf('Seg%d.mat', k)); 
end

if numsegm == 9
 Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg); %,Segs{1,10}.VOCSeg);
elseif numsegm ==8
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg); % Segs{1,9}.VOCSeg); %,Segs{1,10}.VOCSeg);
elseif numsegm == 10
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg);
elseif numsegm == 15
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg);
elseif numsegm == 16
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg);

elseif numsegm == 17
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg);
elseif numsegm == 19
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg);
elseif numsegm == 20
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg);

elseif numsegm == 21
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg);
elseif numsegm == 22
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg,Segs{1,22}.VOCSeg);
elseif numsegm == 23
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg,Segs{1,22}.VOCSeg,Segs{1,23}.VOCSeg);
elseif numsegm == 24
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg,Segs{1,22}.VOCSeg,Segs{1,23}.VOCSeg,Segs{1,24}.VOCSeg);
elseif numsegm == 27
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg,Segs{1,22}.VOCSeg,Segs{1,23}.VOCSeg,Segs{1,24}.VOCSeg,Segs{1,25}.VOCSeg,Segs{1,26}.VOCSeg,Segs{1,27}.VOCSeg);

elseif numsegm == 25
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg,Segs{1,22}.VOCSeg,Segs{1,23}.VOCSeg,Segs{1,24}.VOCSeg,Segs{1,25}.VOCSeg);
elseif numsegm == 29
    Fluxes=vertcat(Segs{1,1}.VOCSeg,Segs{1,2}.VOCSeg,Segs{1,3}.VOCSeg,Segs{1,4}.VOCSeg,Segs{1,5}.VOCSeg,Segs{1,6}.VOCSeg,Segs{1,7}.VOCSeg,Segs{1,8}.VOCSeg, Segs{1,9}.VOCSeg,Segs{1,10}.VOCSeg,Segs{1,11}.VOCSeg,Segs{1,12}.VOCSeg,Segs{1,13}.VOCSeg,Segs{1,14}.VOCSeg,Segs{1,15}.VOCSeg,Segs{1,16}.VOCSeg,Segs{1,17}.VOCSeg,Segs{1,18}.VOCSeg,Segs{1,19}.VOCSeg,Segs{1,20}.VOCSeg,Segs{1,21}.VOCSeg,Segs{1,22}.VOCSeg,Segs{1,23}.VOCSeg,Segs{1,24}.VOCSeg,Segs{1,25}.VOCSeg,Segs{1,26}.VOCSeg,Segs{1,27}.VOCSeg,Segs{1,28}.VOCSeg,Segs{1,29}.VOCSeg);

end
 %Fluxes =(vertcat(Segs{:,:}));

save FluxesLOD_2km Fluxes
writetable(Fluxes,'FluxesLOD_2km.csv');

% output for ArcGIS cannot include NaNs
FluxesGIS = table2array(Fluxes);
FluxesGIS=fillmissing(FluxesGIS,"linear");
FluxesGIS = array2table(FluxesGIS);
FluxesGIS.Properties.VariableNames = Fluxes.Properties.VariableNames;
writetable(FluxesGIS,'FluxesGIS_2km.csv');

%plot(Fluxes.Out_doy,Fluxes.x69_0698776245_C5H9_)
%% LOD
nVOC =(size(Fluxes,2)-15)/6;
FluxLODs=table2array(Fluxes(:,(16+2*nVOC):(15+3*nVOC)));
Fluxtable= table2array(Fluxes(:,16:(15+nVOC)));
LODflagmatrix = abs(Fluxtable) < (2.*FluxLODs); %2 sigma Flux LOD for now
FluxesLODcorr = Fluxtable;
FluxesLODcorr(LODflagmatrix)=0;
FluxesLODcorr =array2table(FluxesLODcorr);
FluxesLODcorr.Properties.VariableNames = Fluxes(:,16:(15+nVOC)).Properties.VariableNames;

%stackedplot(FluxesLODcorr);
%figure; plot(Fluxes.Out_doy,Fluxes.x143_1430358887_C9H19O_,Fluxes.Out_doy,Fluxes.LODx143_1430358887_C9H19O_,Fluxes.Out_doy,2.*Fluxes.LODx143_1430358887_C9H19O_,Fluxes.Out_doy,3.*Fluxes.LODx143_1430358887_C9H19O_); legend('Nonanal','random flux','2x random flux','3x random flux')
% figure; plot(Fluxes.Out_doy,Fluxes.x121_1011734009_C9H13_,Fluxes.Out_doy,Fluxes.LODx121_1011734009_C9H13_,Fluxes.Out_doy,2.*Fluxes.LODx121_1011734009_C9H13_,Fluxes.Out_doy,3.*Fluxes.LODx121_1011734009_C9H13_); legend('TMB','random flux','2x random flux','3x random flux')
% figure; plot(Fluxes.Out_doy,Fluxes.x371_1012268066_C10H31O5Si5_,Fluxes.Out_doy,-Fluxes.LODx371_1012268066_C10H31O5Si5_,Fluxes.Out_doy,-2.*Fluxes.LODx371_1012268066_C10H31O5Si5_,Fluxes.Out_doy,-3.*Fluxes.LODx371_1012268066_C10H31O5Si5_); legend('D5','random flux','2x random flux','3x random flux')


FlagsPerCol = sum(LODflagmatrix,1)./size(LODflagmatrix,1);
LODflagmass = FlagsPerCol == 1;
remainMass = nVOC -sum(LODflagmass);

FluxesLODcorr = horzcat(Fluxes(:,1:15),FluxesLODcorr);

writetable(FluxesLODcorr,'FluxesLODcorr.csv');
%% Plots
% make maps for chosen masses
titles = FluxesLODcorr.Properties.VariableNames;
D = cd;
F = 'MapsFluxLOD';
if ~exist(fullfile(D,F),'dir')
    mkdir(fullfile(D,F))
end
FolderName = (fullfile(D,F));   % Your destination folder

for i=16:(15+nVOC)
    %figure('visible','off');geoscatter(Fluxes.Out_Lat,Fluxes.Out_Lon,20,(table2array(Fluxes(:,i))),'filled')
    fluxnorm = abs(table2array(FluxesLODcorr(:,i))./max(table2array(FluxesLODcorr(:,i))));  
    fluxnorm(fluxnorm==0) = NaN;
    figure('visible','off');geoscatter(FluxesLODcorr.Out_Lat,Fluxes.Out_Lon,80.*fluxnorm,(table2array(Fluxes(:,i))),'filled')
    geobasemap topographic
    %geobasemap satellite
    %geobasemap landcover
    colormap jet
    cb = colorbar();
    
    cb.Label.String = 'mg m^{-2} h^{-1}';
    title(titles{i})
    % caxis([0 scale(i)])
    FigName= sprintf(titles{i});
    saveas(gcf,(fullfile(FolderName,FigName)),'png')  ;
    close all
end

% %% automatically save figures 
% 
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
%     for iFig = 1:length(FigList)
%         FigHandle = FigList(iFig);
%         FigNum = (get(FigHandle, 'Number'));
%       %FigName   = num2str(get(FigHandle, 'Number'));
%       FigName= sprintf(titles{FigNum+13});
%       set(0, 'CurrentFigure', FigHandle);
%       saveas(gcf,(fullfile(FolderName,FigName)),'png')  ;
%     end
% 
% close all