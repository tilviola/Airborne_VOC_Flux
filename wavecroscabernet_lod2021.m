function spec = wavecroscabernet_lod2021(VOC, vw, delay, s_0, s_0_voc, dj, zi, pbl, xkm, Ts,flag)

%change to previous versions: covariance peak is searched for as a local
%max, and covar function is smoothed following Taipale 2010
% this function assumes rotated windspeeds and true vertical velocity, no
% detrending is done
%change Feb 4, 2022: changed dt calculation from mean to mode
%(most common time interval instead of avg time interval), added noise flag
%if cov peak is just noise
% change Feb 23, 2022: added white noise wavelet flux LOD, increased
% search window for cov peak to 4 s

%INPUT
% VOC: 
% column 1: doy
% column 2: VOC (whatever unit)

% vw:
% column 1: time of day in hours
% column 2: vertical windspeed in m/s
% column 3: temperature

% delay time between VOC dataset and vertical velocity in seconds: if VOC
% timestamp lags behind delay time should be positive
% s_0: 2 x sampling rate of 10 Hz sensor
% s_0_voc: 2 x sampling rate of disjunct VOC sensor
% dj: houskeeping values for wavelet analysis - use 0.25 as default

% zi: measurement height in km (needed for systematic error)
% pbl: boundary layer height in km (needed for random error)
% xkm: aircraft velocity in km/s
% Ts: integral time scale in s (used for error calculation)
% flag: determines if contour plot is made or not

% OUTPUT
% spec structure: global spectra (T, w, and wT) etc as following
%  spec.fqfft: fft frequency (s-1) for w, T and w'T'
%  spec.specTfft: fft T spectrum in f x Spec
%  spec.specwfft: fft w spectrum in f x Spec
%  spec.wT: fft w'T' spectrum in f x Spec
%  spec.scalewT: wavelet scale
%  spec.wavewT: (variancew*varianceT)*sum(wc_wT')/nw;
%  spec.fqwave=1./periodw; wavelet frequency
%  spec.wave=wc_wT; global wavelet wT spectrum
%  spec.scale=scalew; wavelet w scale (normally the same as ScaleT and
%  scale wT
%  spec.coi: cone of influence for highfrequency data
%  spec.coivoc: cone of incfluence for lowfrequency (disjunct) data
%  spec.timx: VOC time stamp
%  spec.wavefluxvoc: VOC waveletflux
%  spec.fftfluxvoc: VOC fft flux
%  spec.se: systematic error
%  spec.re: random error
%  spec.de: error due to disjunct sampling
% spec.wavelet_noise_mean: mean random flux at +/-175 to 180 s lag time
% spec.wavelet_noise_std: standard deviation of flux at +/-175 to 180 s lag time
% spec.fluxSNRwave: signal to noise ratio of wavelet flux

%vwi=interp1(vw(:,1),vw(:,2),VOC(:,1));
%VOCi=interp1(VOC(:,1)*3600-delay,VOC(:,2),vw(:,1)*3600);
VOCi=interp1(VOC(:,1)*3600*24-delay,VOC(:,2),vw(:,1)*3600*24); %VOC conc, time corrected for delay
indx=(~isnan(VOCi));
%dt = mode(diff(VOC(:,1)*3600*24)) %time interval of disjunct VOC sensor

%----FFT covariance of vw and VOC----------
[vocflux,lagc]=xcov(vw(indx,2),VOCi(indx),3000); %covariance of vw and VOC searching in a window of 3000 sampling intervals
center = round(length(vocflux)./2); %define center to find local maximum of covariance
vocfluxf=filtfilt(ones(10,1)/10,1,vocflux); % smoothing the covariance peak so that it doesn't get influenced by noise in itself
% this is following Taipale et al., 2010

centerwindow =40;
if abs(min(vocfluxf((center-centerwindow):(center + centerwindow)))) > abs(max(vocfluxf((center-centerwindow):(center + centerwindow))))
    [fftflux,lagtime]=min(vocfluxf((center-centerwindow):(center + centerwindow))); %returns min, ignoring Nans -> maximum covariance divided by no. of points is flux
     fftflux= fftflux/size(VOCi,1);
     lagtime=(lagtime-(centerwindow+1))/10;
    vocfluxn=vocfluxf./fftflux;
else [fftflux,lagtime] = max(vocfluxf((center-centerwindow):(center + centerwindow))); %returns max, ignoring Nans -> maximum covariance divided by no. of points is flux
    fftflux= fftflux/size(VOCi,1);
    lagtime=(lagtime-(centerwindow+1))/10;
    vocfluxn=vocfluxf./fftflux;
end
%figure; plot(lagc,vocflux); %plots the crossvariance peak of vw and VOC vs lagtime
lagc1=lagc/10; % in s
figure; plot(lagc1-lagtime,vocfluxn); xlabel("lag time (s)"); ylabel("w'c'")

%calculate flux errors
spec.se=2.2*pbl*sqrt(zi/pbl)/(xkm*(vw(size(vw,1),1)-vw(1,1))*3600*24);
spec.re=1.75*(zi/pbl)^(0.25)*sqrt(pbl/(xkm*(vw(size(vw,1),1)-vw(1,1))*3600*24));

% disjunct sampling
dt_dis=mode(diff(VOC(:,1)))*3600*24;
Tsample=(VOC(size(VOC,1),1)-VOC(1,1))*3600*24;
spec.de=dt_dis/Tsample*((cosh(dt_dis/Ts/2)/sinh(dt_dis/Ts/2))-dt_dis/Tsample*(1-exp(-Tsample/Ts))/(2*sinh(dt_dis/Ts/2)^2));

if 1
 %global wT fft spectrum
    
  %[wTfft,lagwT]=xcov(vw(:,2),vw(:,3));

fftvw=fft(vw(:,2));
fftvw2=fft(vw(indx,2));
fftT=fft(vw(:,3));



dt10hz=mode(diff(vw(:,1)*3600*24));
step=abs(vw(size(vw,1),1)-vw(1,1))*3600*24/dt10hz;
fq=(((1:floor((step)/2)))/step)./dt10hz; %floor: rounding to next lower integer, what is step doing=? 
vw2(:,1)=vw(indx,1);
vw2(:,2)=vw(indx,2);
step2=abs(vw2(size(vw2,1),1)-vw2(1,1))*3600*24/dt10hz;
fq2=(((1:floor((step2)/2)))/step2)./dt10hz;

fftvoc=fft(VOCi(indx));
fftwvoc=fft(vw2(:,2));

specT=real(fftT(2:ceil(size(fq,2)+1))).^2+imag(fftT(2:ceil(size(fq,2)+1))).^2;
specw=real(fftvw(2:ceil(size(fq,2)+1))).^2+imag(fftvw(2:ceil(size(fq,2)+1))).^2;

specvoc=real(fftvoc(2:ceil(size(fq2,2)+1))).^2+imag(fftvoc(2:ceil(size(fq2,2)+1))).^2;
specT=specT/step/2;
specw=specw/step/2;
specvoc=specvoc/step2/2;

wT=(real(fftT(2:ceil(size(fq,2)+1))).*real(fftvw(2:ceil(size(fq,2)+1)))+imag(fftT(2:ceil(size(fq,2)+1))).*imag(fftvw(2:ceil(size(fq,2)+1))))/step/2;
wVOC=(real(fftvoc(2:ceil(size(fq2,2)+1))).*real(fftwvoc(2:ceil(size(fq2,2)+1)))+imag(fftvoc(2:ceil(size(fq2,2)+1))).*imag(fftwvoc(2:ceil(size(fq2,2)+1))))/step2/2;

nn=15; %bin spectra
[f_w,y_w] = binLogSpaced(fq,fq'.*specw,nn);
[f_T,y_T] = binLogSpaced(fq,fq'.*specT,nn);
[f_wT,y_wT] = binLogSpaced(fq,fq'.*wT,nn);
[f_wvoc,y_wvoc] = binLogSpaced(fq2,fq2'.*wVOC,nn);

% flux detection limit: flux noise criterium using STD noise of covariance between +/- 160-180s (Spirig et al. 2005)
                        idx = [(-240/0.1):(-220/0.1), (220/0.1):(240/0.1)];
                        flux_noise_mean = abs(mean(vocflux(center+idx)/size(VOCi,1),'omitnan')); 
                        flux_noise_std = std(vocflux(center+idx)/size(VOCi,1),'omitnan');

                        if abs(3*flux_noise_std) >= abs(fftflux)
                            lagtime=0; %if the covariance peak is just noise, use standard lag time!!
                            noiseflag=1;
                        else
                            noiseflag =0;
                        end

% flux detection limit: flux noise criterium using RMSE noise of covariance
% between +/- 160-180s (Langford 2015)
                        idx_left = (-240/0.1):(-220/0.1);
                        idx_right = (220/0.1):(240/0.1);
                        flux_noise_rmse = sqrt( 0.5*( std(vocflux(center+idx_left),'omitnan')^2 + mean(vocflux(center+idx_left),'omitnan')^2 + std(vocflux(center+idx_right),'omitnan')^2 + mean(vocflux(center+idx_right),'omitnan')^2 ) )/size(VOCi,1);

% flux detection limit: random shuffle criterium as described by Billesbach
% 2011 - Langford 2015 says this overestimates the uncertainty
                        for irand = 1:10
                            w_rand = vw((randperm(length(vw(indx,2)))),2);
                            xcov_rand(irand) = xcov(w_rand, VOCi(indx), 0, 'unbiased');
                        end
                        random_flux = std(xcov_rand,'omitnan');

  % estimate white noise using autocovariance (Lenschow 2000, Mauder 2013, Langford 2015)
                        autocov_c = xcov(VOCi(indx), VOCi(indx), 10, 'unbiased');
                        white_noise_c = ( sqrt( autocov_c(11) ) - sqrt( polyval( polyfit( 1:9, (autocov_c(1:9))', 1 ), 11 ) ) );
                        random_error_noise = white_noise_c * std(vw(indx,2),'omitnan') / sqrt(size(VOCi,1));

%%---FFT wT flux ---
[wTflux,~]=xcov(vw(:,2),vw(:,3),4000);

center = ceil(length(wTflux)./2); %define center to find local maximum of covariance
wTflux=filtfilt(ones(20,1)/20,1,wTflux); %smooth so that the maximum is not impeded by noise

if abs(min(wTflux((center-20):(center + 20)))) > abs(max(wTflux((center-20):(center + 20))))
    [wTfftflux,lagtimeT]=min(wTflux((center-20):(center + 20)));
    wTfftflux=wTfftflux./size(vw,1); % maximum covariance divided by no. of points is flux
else [wTfftflux,lagtimeT]=max(wTflux((center-20):(center + 20))); %returns max, ignoring Nans -> maximum covariance divided by no. of points is flux
    wTfftflux=wTfftflux./size(vw,1);
end



spec.fqfft=f_w;
spec.specTfft=y_T;
spec.specwfft=y_w;
spec.wT=y_wT;
spec.wvoc=y_wvoc;
spec.fftfluxwT=wTfftflux; %fft flux for wT
spec.flux_noise_std = flux_noise_std;
spec.flux_noise_mean = flux_noise_mean;
spec.fluxSNRfft = (abs(fftflux) - abs(spec.flux_noise_mean))./abs(spec.flux_noise_std);
spec.flux_noise_rmse =flux_noise_rmse  ;
spec.random_flux =random_flux;
spec.random_error_noise =random_error_noise;
spec.lagtime=lagtime;
end

%% wT wavelet
sstw=vw(:,2);
sstT=vw(:,3);

variancew = std(sstw)^2;
sstw = (sstw - mean(sstw))/sqrt(variancew) ;
varianceT = std(sstT)^2;
sstT = (sstT - mean(sstT))/sqrt(varianceT) ;

nw = length(sstw);
timew = [0:length(sstw)-1]*dt10hz + 0 ;  % construct time array
xlimw = [0,max(timew)*xkm];  % plotting range
padw = 1;      % pad the time series with zeroes (recommended)
dj;% = 0.2;    % this will do x sub-octaves per octave
s0 = s_0;    % 2 x sampling interval
j1 = round(log2(round(size(sstw,1)/ceil(32/(s_0_voc/s_0)))))/dj;    % e.g. use log2 of half the sample size with dj sub-octaves each - scale 64 by VOC time interval and 10 Hz 
lag1 = 0.72;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wavew,periodw,scalew,coiw] = WAVELET(sstw,dt10hz,padw,dj,s0,j1,mother);
[waveT,periodT,scaleT,coiT] = WAVELET(sstT,dt10hz,padw,dj,s0,j1,mother);

%power = (abs(wave)).^2 ;        % compute wavelet power spectrum

%crossspectrum
wc_wT=real(wavew).*real(waveT)+imag(wavew).*imag(waveT);

% Significance levels: (variance=1 for the normalized SST)
[signif_wT,fft_theor_wT] = wave_signif(1.0,dt10hz,scalew,0,lag1,-1,-1,mother);
sig95_wT = (signif_wT')*(ones(1,nw));  % expand signif --> (J+1)x(N) array
sig95_wT = wc_wT./ sig95_wT;         % where ratio > 1, power is significant


avg = find((scalew >= min(scalew)	& (max(scalew))));%	& (scale69 <= 3600)));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scalew')*(ones(1,nw));  % expand scale --> (J+1)x(N) array
scale_avg = wc_wT./ scale_avg;   % [Eqn(24)]
scale_avg = sqrt(variancew*varianceT)*dj*dt10hz/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
scaleavg_signif = wave_signif(sqrt(variancew*varianceT),dt10hz,scalew,2,lag1,-1,[min(scalew),max(scalew)],mother);

spec.wavefluxwT=scale_avg'; % write wavelet wT flux to structure
clear avg 
clear scale_avg
clear scaleavg_signif

spec.scalewT=scalew;
spec.wavewT=(variancew*varianceT)*sum(wc_wT,2)/nw;
%spec.wavewTub=(variancew*varianceT)*sum(wc_wT_ub,2)/nw;
spec.fqwave=1./periodw;
spec.wave=wc_wT;
spec.coi=coiw;
spec.timew=timew; % time stamp for high frequency data

avg = (scalew >= min(scalew)	& (max(scalew)));%	& (scale69 <= 3600)));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scalew')*(ones(1,nw));  % expand scale --> (J+1)x(N) array
scale_avg = wc_wT./ scale_avg;   % [Eqn(24)]
scale_avg = sqrt(variancew*varianceT)*dj*dt10hz/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
scaleavg_signif = wave_signif(sqrt(variancew*varianceT),dt10hz,scalew,2,lag1,-1,[min(scalew),max(scalew)],mother);
wTwavelet=scale_avg';
spec.wT_time=wTwavelet;
spec.wstar=(9.81*pbl*1000*abs(wTwavelet./vw(:,3))).^(1/3);
%plot wT wavelet spectrum
if flag
    figure(3)
%--- Plot time series
subplot('position',[0.1 0.75 0.65 0.2])
plot(timew*xkm,sstw,'b', timew*xkm,sstT, 'r')
set(gca,'XLim',xlimw(:))
ylabel('[c - mean(c)]/sqrt(c)')
title('a) w, T')
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(periodw))):fix(log2(max(periodw))));
%sig95(find(sig95<=0))=1e-5;
surface(timew*xkm,log2(periodw),wc_wT/max(max(abs(wc_wT)))); shading interp; colorbar
%contourf(time,log2(period59),log2(wc),log2(levels));  %*** or use 'contourfill'
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
%xlabel('Time (s)')
ylabel('Period (s)')
title('b) Wavelet Cross Spectrum')
set(gca,'XLim',xlimw(:))
set(gca,'YLim',log2([min(periodw),max(periodw)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(timew*xkm,log2(periodw),sig95_wT,[-99,1],'k'); % or km? speed x time
colormap((jet))
hold on
% cone-of-influence, anything "below" is dubious
plot3(timew*xkm,log2(coiw),ones(1,size(coiw,2)).*max(max(wc_wT/max(max(abs(wc_wT))))),'-.w')
hold off

%--- Plot scale-average time series
subplot('position',[0.1 0.07 0.65 0.2])
halfhour(1:size(timew,2))=mean(wTwavelet);
%ghh=max(cumsum(global_ws(1:32)))/32;
plot(timew*xkm,wTwavelet,'r');
set(gca,'XLim',xlimw(:))
xlabel('Distance [km]')
ylabel('Av. crossvariance ')
title('c) Scale-average Time Series')
hold on
%plot(time,halfhour,'--')
%plot(xlim,scaleavg_signif+[0,0],'--')
hold off

end

clear scale_avg
clear scaleavg_signif

%% VOC wavelet
clear indx

vwi=interp1(vw(:,1)*3600*24+delay-lagtime,vw(:,2),VOC(:,1)*3600*24);


indx=(~isnan(vwi));

VOCc = VOC(indx,2); 
vwic = vwi(indx); 

%------------------------------------------------------ Computation

% normalize by standard deviation 
variance1 = std(VOCc)^2;
VOCc = (VOCc - mean(VOCc))/sqrt(variance1) ;
variance2 = std(vwic)^2;
vwic = (vwic - mean(vwic))/sqrt(variance2) ;

n = length(VOCc);
dt = mode(diff(VOC(:,1)*3600*24));
time = [0:length(VOCc)-1]*dt + 0 ;  % construct time array
xlim = [0,max(time)*xkm];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
s0 = s_0_voc;    % e.g. s_0=24 this says start at a scale of x seconds
j1 = round(log2(round(size(VOCc,1))/32))/dj;    % e.g. use log2 of half the sample size with dj sub-octaves each 64 = scaling factor to choose the right wavelet scales
lag1 = 0.72;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
%[wave69,period69,scale69,coi69] = WAVELET(sst69,dt,pad,dj,s0,j1,mother);
[wave69,period69,scale69,coi69] = WAVELET(VOCc,dt,pad,dj,s0,j1,mother);
[wave71,period71,scale71,coi71] = WAVELET(vwic,dt,pad,dj,s0,j1,mother);
%power = (abs(wave)).^2 ;        % compute wavelet power spectrum
spec.coivoc=coi69;
spec.scalevoc=scale69;

%crossspectrum
wc6971=real(wave69).*real(wave71)+imag(wave69).*imag(wave71);
spec.wavewVOC=(variance1*variance2)*sum(wc6971,2)/n;
spec.fqwavevoc=1./period69;
spec.crosspecVOC = wc6971;

% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale69,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = wc6971./ sig95;         % where ratio > 1, power is significant
%figure; surface(sig95); shading interp

            if 0
            % Global wavelet spectrum & significance levels:
            global_ws = sqrt(variance1*variance2)*(sum(wc6971')/n);   % time-average over all times
            dof = n - scale69;  % the -scale corrects for padding at edges
            global_signif = wave_signif(sqrt(variance2*variance1),dt,scale69,1,lag1,-1,dof,mother);
            
            % Scale-average between El Nino periods of 2--8 years
            
            avg = find((scale69 >= 1	& (scale69 <= 30)));
            Cdelta = 0.776;   % this is for the MORLET wavelet
            scale_avg = (scale69')*(ones(1,n));  % expand scale --> (J+1)x(N) array
            scale_avg = wc6971./ scale_avg;   % [Eqn(24)]
            scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
            scaleavg_signif = wave_signif(sqrt(variance1*variance2),dt,scale69,2,lag1,-1,[0.12,16],mother);
            
            sf(:,1)=scale_avg';
            
            avg = find((scale69 >= max(scale69)/2	& (max(scale69))));
            Cdelta = 0.776;   % this is for the MORLET wavelet
            scale_avg = (scale69')*(ones(1,n));  % expand scale --> (J+1)x(N) array
            scale_avg = wc6971./ scale_avg;   % [Eqn(24)]
            scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
            scaleavg_signif = wave_signif(sqrt(variance1*variance2),dt,scale69,2,lag1,-1,[ max(scale69)/2,(max(scale69))],mother);
            
            sf(:,2)=scale_avg';
            end

Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scale69')*(ones(1,n));  % expand scale --> (J+1)x(N) array
scale_avg = wc6971./ scale_avg;   % [Eqn(24)]
scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg((scale69 >= min(scale69)	& (max(scale69))),:));   % [Eqn(24)]
scaleavg_signif = wave_signif(sqrt(variance1*variance2),dt,scale69,2,lag1,-1,[min(scale69),max(scale69)],mother);
sf(:,3)=scale_avg';

%output wavelet
waveflux = scale_avg';
spec.timx=time; %VOC time stamp
spec.wavefluxvoc=waveflux; %VOC waveletflux
spec.fftfluxvoc=fftflux; % VOC fft flux
spec.scaleavg_signif =scaleavg_signif;

  %% COI influence
period_big = repmat(period69',1,length(time));
coi_big    = repmat(coi69,length(scale69),1);
icoi = (coi_big>period_big);
crosspeccoi = (wc6971.*~icoi);
ratiocoi = sum(crosspeccoi)./sum(wc6971);
coi_avg = (scale69')*(ones(1,n));
coi_avg = crosspeccoi./ coi_avg;   % [Eqn(24)]
coi_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(coi_avg);   % [Eqn(24)]
qcoi = abs(coi_avg./scale_avg) ; %fraction of flux under COI
spec.qcoi = qcoi;
spec.qcoir = ratiocoi;

%% Waveflux LOD

%flux detection limit: flux noise criterium using STD noise of covariance between +/- 160-180s (Spirig et al. 2005)
                        delayLOD = [(-180):(-175),(175):(180)]'; %longer because large eddies. for sake of calculation time, use less

                        for i = 1:length(delayLOD)
                            LOD=[];
                                            
                        vwiLOD=interp1((vw(:,1)*3600*24+delayLOD(i)-lagtime),vw(:,2),VOC(:,1)*3600*24);
                        LOD.w=vwiLOD(~isnan(vwiLOD));
                        LOD.x = VOC(~isnan(vwiLOD),2);
                        LOD.t = VOC(~isnan(vwiLOD),1)*3600*24+delayLOD(i)-lagtime;

                        dt1 = mode(diff(LOD.t));
                        j11 = round(log2(round(length(LOD.w))/32))/dj;    
                        n = length(LOD.t);

                        % normalize by standard deviation 
                        variance11 = std(LOD.x)^2;
                        LOD.x = (LOD.x - mean(LOD.x))/sqrt(variance11) ;
                        variance21 = std(LOD.w)^2;
                        LOD.w = (LOD.w - mean(LOD.w))/sqrt(variance21) ;

                            if n > 1800
                            %Wavelet transformation
                            [wave691,period691,scale691,coi691] = WAVELET(LOD.x,dt1,1,dj,s_0_voc,j11,mother);
                            [wave711,period711,scale711,coi711] = WAVELET(LOD.w,dt1,1,dj,s_0_voc,j11,mother);
    
                          
                            %crossspectrum
                            wc69711=real(wave691).*real(wave711)+imag(wave691).*imag(wave711);
                                                 
                                                    
                            % Scale-average 
                            
                            Cdelta = 0.776;   % this is for the MORLET wavelet
                            scale_avg = (scale691')*(ones(1,n));  % expand scale --> (J+1)x(N) array
                            scale_avg = wc69711./ scale_avg;   % [Eqn(24)]
                            scale_avg = sqrt(variance11*variance21)*dj*dt/Cdelta*sum(scale_avg((scale691 >= min(scale691)	& (max(scale691))),:));   % [Eqn(24)]
    
                            %output wavelet
                            LODflux{i} = scale_avg';
                            LODfluxes=vertcat(LODflux{:});
   
                          spec.wavelet_noise_mean = mean(abs(LODfluxes),'all','omitnan'); 
                          spec.wavelet_noise_std = abs(std(LODfluxes,0,'all','omitnan')); 
                          if dt <= 0.15 %for 10 Hz
                          spec.wavelet_noise_movstd = abs(std(movmean(LODfluxes,400,'omitnan','Endpoints','discard'))); 
                          elseif dt > 0.15  %for 2 Hz
                              spec.wavelet_noise_movstd = abs(std(movmean(LODfluxes,80,'omitnan','Endpoints','discard'))); 
                          end
                            else 
                              spec.wavelet_noise_mean = 3*flux_noise_mean; 
                              spec.wavelet_noise_std = 3*flux_noise_std; 
                              spec.wavelet_noise_movstd = 3*flux_noise_std;
                            end
                        end
                        
      spec.fluxSNRwave = abs((waveflux) - spec.wavelet_noise_mean)./spec.wavelet_noise_std;
spec.noiseflag = noiseflag;

%% Waveflux LOD with white noise
% generate white noise with same stdev and mean and sample number as the
% VOC timeseries
LODw.x = VOC(~isnan(vwi),2);
LODw.xsdv = std(LODw.x);
LODw.xmean = mean(LODw.x);
%autocov_c = xcov(VOC(~isnan(vwi),2));

LODw.w = vwi(~isnan(vwi));
LODw.wsdv = std(LODw.w);
LODw.wmean = mean(LODw.w);
LODw.t = VOC(~isnan(vwi),1)*3600*24-lagtime;
n = length(LODw.t);
dt1 = mode(diff(LODw.t));
j11 = round(log2(round(length(LODw.w))/32))/dj;

for irand = 1:10
    LODw.xwhite = white_noise_c.*randn(n,1) ;

    % normalize by standard deviation
    variance11 = std(LODw.xwhite)^2;
    LODw.xwhite = (LODw.xwhite - mean(LODw.xwhite))/sqrt(variance11) ;
    variance21 = std(LODw.w)^2;
    LODw.w = (LODw.w - mean(LODw.w))/sqrt(variance21) ;


    %Wavelet transformation
    [wave691,period691,scale691,coi691] = WAVELET(LODw.xwhite,dt1,1,dj,s_0_voc,j11,mother);
    [wave711,period711,scale711,coi711] = WAVELET(LODw.w,dt1,1,dj,s_0_voc,j11,mother);


    %crossspectrum
    wc69711=real(wave691).*real(wave711)+imag(wave691).*imag(wave711);


    % Scale-average

    Cdelta = 0.776;   % this is for the MORLET wavelet
    scale_avg = (scale691')*(ones(1,n));  % expand scale --> (J+1)x(N) array
    scale_avg = wc69711./ scale_avg;   % [Eqn(24)]
    scale_avg = sqrt(variance11*variance21)*dj*dt/Cdelta*sum(scale_avg((scale691 >= min(scale691)	& (max(scale691))),:));   % [Eqn(24)]

    %output wavelet
    LODfluxw{irand} = scale_avg';
end
                            
                      
                                  
                LODfluxesw=vertcat(LODfluxw{:});
                spec.wavelet_white_mean = mean(abs(LODfluxesw),'all','omitnan'); 
                spec.wavelet_white_std = abs(std(LODfluxesw,0,'all','omitnan'));

                if dt <= 0.15 %for 10 Hz
                    spec.wavelet_white_movstd = abs(std(movmean(LODfluxesw,400,'omitnan','Endpoints','discard')));
                elseif dt > 0.15  %for 2 Hz
                    spec.wavelet_white_movstd = abs(std(movmean(LODfluxesw,80,'omitnan','Endpoints','discard')));
                end

                % the random flux cannot be smaller than the random FFT
                % flux, so finally using the maximum between Spirig and this white
                % noise approach
spec.waveflux_sigma=max([spec.wavelet_white_movstd;spec.flux_noise_std]);

%% Waveflux LOD with random shuffle
% flux detection limit: random shuffle criterium as described by Billesbach 2011
% according to Langford 2015, this overestimates the LOD. I can confirm
% that.
%                         for irand = 1:10
%                             w_rand = vwi((randperm(length(vwi))));
%                             LODs=[];
%                                             
%                         
%                         LODs.w=w_rand(~isnan(w_rand));
%                         LODs.x = VOC(~isnan(w_rand),2);
%                         LODs.t = VOC(~isnan(w_rand),1)*3600*24-lagtime;
% 
%                         dt1 = mode(diff(LODs.t));
%                         j11 = round(log2(round(length(LODs.w))/32))/dj;    
%                         n = length(LODs.t);
% 
%                         % normalize by standard deviation 
%                         variance11 = std(LODs.x)^2;
%                         LODs.x = (LODs.x - mean(LODs.x))/sqrt(variance11) ;
%                         variance21 = std(LODs.w)^2;
%                         LODs.w = (LODs.w - mean(LODs.w))/sqrt(variance21) ;
% 
%                             
%                             %Wavelet transformation
%                             [wave691,period691,scale691,coi691] = WAVELET(LODs.x,dt1,1,dj,s_0_voc,j11,mother);
%                             [wave711,period711,scale711,coi711] = WAVELET(LODs.w,dt1,1,dj,s_0_voc,j11,mother);
%     
%                           
%                             %crossspectrum
%                             wc69711=real(wave691).*real(wave711)+imag(wave691).*imag(wave711);
%                                                  
%                                                     
%                             % Scale-average 
%                             
%                             Cdelta = 0.776;   % this is for the MORLET wavelet
%                             scale_avg = (scale691')*(ones(1,n));  % expand scale --> (J+1)x(N) array
%                             scale_avg = wc69711./ scale_avg;   % [Eqn(24)]
%                             scale_avg = sqrt(variance11*variance21)*dj*dt/Cdelta*sum(scale_avg((scale691 >= min(scale691)	& (max(scale691))),:));   % [Eqn(24)]
%     
%                             %output wavelet
%                             LODfluxs{irand} = scale_avg';
%                 
%                             
%                         end
%                                   
%                 LODfluxess=vertcat(LODfluxs{:});
%                 spec.wavelet_shuff_mean = mean(abs(LODfluxess),'all','omitnan'); 
%                 spec.wavelet_shuff_std = abs(std(LODfluxess,0,'all','omitnan')); 
%% plotting
if flag
%------------------------------------------------------ Plotting
figure(4)
%--- Plot time series
subplot('position',[0.1 0.75 0.65 0.2])
plot(time*xkm,VOCc,'b', time*xkm,vwic, 'r')
set(gca,'XLim',xlim(:))
ylabel('[c - mean(c)]/sqrt(c)')
title('a) VOC')
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period69))):fix(log2(max(period69))));
%sig95(find(sig95<=0))=1e-5;
surface(time*xkm,log2(period69),wc6971/max(max(abs(wc6971)))); shading interp; colorbar
%contourf(time,log2(period59),log2(wc),log2(levels));  %*** or use 'contourfill'
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
%xlabel('Time (s)')
ylabel('Period (s)')
title('b) Wavelet Cross Spectrum')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period69),max(period69)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(time*xkm,log2(period69),sig95,[-99,1],'k'); % or km? speed x time
colormap((jet))
hold on
% cone-of-influence, anything "below" is dubious
plot3(time*xkm,log2(coi69),ones(1,size(coi69,2)).*max(max(wc_wT/max(max(abs(wc_wT))))),'-.w')
hold off

%--- Plot scale-average time series
subplot('position',[0.1 0.07 0.65 0.2])
halfhour(1:size(time,2))=mean(sf(:,1));
%ghh=max(cumsum(global_ws(1:32)))/32;
plot(time*xkm,sf(:,1),'r',time*xkm,sf(:,2),'b',time*xkm,sf(:,3),'g');
set(gca,'XLim',xlim(:))
xlabel('Distance [km]')
ylabel('Av. crossvariance')
title('c) Scale-average Time Series')

hold on
%plot(time,halfhour,'--')
%plot(xlim,scaleavg_signif+[0,0],'--')
hold off

end



      
 
  