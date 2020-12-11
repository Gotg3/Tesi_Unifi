
  addpath(genpath('C:\Users\Gianluca\Desktop\MATLAB\SIMULATION_01'));
  addpath('C:\Users\Gianluca\Desktop\MATLAB\SIMUL_SBB\SBB_PULSED');

  run('C:\Users\Gianluca\Desktop\MATLAB\SIMULATION_01\Lettura file\leggiSBB.m');

clear vel

leggiSBB;


%%%%%%%%%%%%%%%%%%%%%%%%RICERCA DEI PARAMETRI NEL NOME%%%%%%%%%%%%%%%%%%%%%

%filename è in leggiSBB
prfIndex=strfind(filename, 'PRF');
FtxIndex=strfind(filename, 'f');  %ricerca indice della lettera f
Sxzindex=strfind(filename, 'xz');    %ricerca indice lettera z
Syzindex=strfind(filename,'y');
prf=str2double(filename(prfIndex(1)+3:prfIndex(1)+4));
PRF=prf*1e2;
Steeringxz=str2double(filename(Sxzindex(1)+2:Sxzindex(1)+3));  %converte in numero il carattere all'idice dopo zeta
Steeringyz=str2double(filename(Syzindex(1)+2:Syzindex(1)+3));
%considero due caratteri perchè il numero è con due cifre, si mette +2
%perchè ha preso l'indice di x
f=str2double(filename(FtxIndex(1)+1));
Ftx=f*1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L,NFFT]=size(sgn);         %dimensioni sgn
c=1480;             %velocità del suono

dGate=c/(2*fs);     %lunghezza gate

firstDp=(c*tt(1)/2);            %prima depth
dp_axis=abs(((0:L-1)*dGate)+firstDp);   %creazione asse
z=dp_axis/sqrt(1+((tand(Steeringxz))^2+(tand(Steeringyz))^2));  %conversione in z


NFFTpart=256;               %nfft parziale

Overlap=0.95;
threshold=-6;               %soglia espressa in db
thPlot=-35;                 %soglia plot

%%%%%%%%%%%%%%%%%%%%%ROUTINE OVERLAP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  NFFToverlap=ceil(Overlap*NFFTpart);         %ceil arrotonda all'intero più vicino (>=)

distOverlap=NFFTpart-NFFToverlap;              %distanza tra una fft e l'altra

%prima limitazione overflow

numberFFT=length(0:distOverlap:NFFT-NFFTpart);% conta il +1 anche con lo 0

FFTmat=zeros(L,NFFTpart);     %matrice di appoggio per memorizzare i valori della fft
firstNFFT=1;                  %indice iniziale fft overlap
lastNFFT=firstNFFT+NFFTpart-1;    %indice finale fft overlap
overflow=ceil(NFFTpart+(numberFFT-1)*distOverlap-NFFT);               %mi dice di quanto esco dai miei NFFT

if overflow>(NFFTpart/2)
    numberFFT=numberFFT-1;     %ammetto solo un overflow di NFFTpart/2
end
Newrow=ones(L,overflow);

fd_axis=(0:NFFTpart-1)/(NFFTpart)*PRF;            %generazione asse frequenze
fd_axis=fd_axis-PRF/2;

FFTcube=zeros(L,NFFTpart,numberFFT);              %matrice di storage delle fft
SumFFT=zeros(L,NFFTpart);

% %%APPLICATA MODIFICA MIA E DI CLAUDIO
% tSpt=zeros(1,numberFFT+2);  % forse qua si doveva aggiungere un ultimo
% % indice per quel discorso della linea, non entra la matrice ma il suo
% % inizio si e probabilmente anche l'ultimo elemento di sgn, è per questo
% % che c'è +2 probabilmente
 tSpt=zeros(1,numberFFT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:numberFFT
    tSpt(k)=firstNFFT;
    %   selezione una parte delle righe con (n:m,:) selezione da n a m
    MatIQ=sgn(:,firstNFFT:lastNFFT);
    FFTcube(:,:,k)=fft(MatIQ,[],2);
    firstNFFT=firstNFFT+distOverlap;
    lastNFFT=firstNFFT+NFFTpart-1;   %aggiornamento indici nfft
end
%modifiche fatte con claudio, corrispondono a quello che avevamo detto,

tSpt(end-1)=tSpt(end-2)+NFFTpart;
tSpt(end)=size(sgn,2);
FFTcube=FFTcube/max(abs(FFTcube(:)));
FFTmed=mean(abs(FFTcube),3);
FFTmed=fftshift(FFTmed,2);

m_fft=20*log10(FFTmed);
norm_fft=m_fft-max(m_fft(:));

%%% plot Multigate Spectral Doppler medio
h=figure(201);
a=axes('parent',h);
imagesc(a,fd_axis,z*1e3,norm_fft,[thPlot 0])
xlabel(a,'Doppler frequency [Hz]')
ylabel(a,'Depth [mm]')
colormap(a,hot)
%%


[~, vesC]=min(abs(z-20e-3));            %ricerca centro vaso

%%%%%%%%%%%%%%%%%%%%%%ROUTINE PARAMETRI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linThreshold=10^(threshold/20);

norm_shift=FFTmed(vesC,:);   %assegna il valore della fft nel centro

% [limits]=findLimitsV2(norm_shift/max(norm_shift), linThreshold);         %uso la function per trovare i limiti
norm_shift_smoothed=10.^(smooth(log10(norm_shift/max(norm_shift))));
[limits]=findLimitsV2(norm_shift_smoothed/max(norm_shift_smoothed), linThreshold);         %uso la function per trovare i limiti

band=abs(diff(fd_axis(limits)));
avrgFreq=sum(fd_axis(limits(1):limits(2)).*norm_shift(limits(1):limits(2))/sum(norm_shift(limits(1):limits(2))));



%%%%%%%%%%SPETTROGRAMMA%%%%%%%%%%%%%
h=figure(200);
a=axes('parent',h);
spt=(squeeze(FFTcube(vesC,:,:)));
spt=fftshift(spt/max(abs(spt(:))),1);
tSpt=tSpt/PRF;
imagesc(a,tSpt*1e3,fd_axis,20*log10(abs(spt)),[thPlot,0]);
hold(a,'all');
plot(a,[tSpt(1),tSpt(end)]*1e3,[0,0],'w','linewidth',2);
hold(a,'off');
set(a, 'YDir', 'normal')
xlabel(a,'Time [ms]');
ylabel(a,'Doppler frequency [Hz]');
colormap(a,hot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vel=avrgFreq*c/(2*Ftx*cosd(90+Steeringxz)); %calcolo della velocità

y=(20*log10(norm_shift));
y=y-max(y);
ysmoothed=20*log10(norm_shift_smoothed);

%%%%Plot
h=figure(101);
a=axes('parent',h);
plot(a,fd_axis,y)
hold(a,'all');
plot(a,fd_axis,ysmoothed)
ylim(a,[thPlot,0]);
xlim(a,[min(fd_axis),max(fd_axis)]);
scatter(a,fd_axis([limits(1),limits(2)]),ysmoothed([limits(1),limits(2)]),'filled');
plot(a,[avrgFreq,avrgFreq],get(gca,'ylim'),'--k');
hold(a,'off');
xlabel(a,'Doppler frequency [Hz]');
ylabel(a,'Magnitude [dB]');

%%%% Risultati
fprintf('------\n Band: %.1f Hz\n Fd: %.1f Hz\n Flow speed: %.1f cm/s\n------\n',band,avrgFreq,vel*100);


%%%%Video
h=figure(1000);
FR=1/(tSpt(2)-tSpt(1));
FRvideo=min(60,FR);
[~,videoName,~]=fileparts(filename);
writerObj=VideoWriter([videoName,'.avi'],'MPEG-4');
writerObj.FrameRate=FRvideo;
writerObj.Quality=85;
open(writerObj);
a=axes('parent',h);
colormap(a,hot);
for i=round(1:FR/FRvideo:size(FFTcube,3))
    imagesc(a,fd_axis,z*1e3,fftshift(20*log10(abs(FFTcube(:,:,i))),2),[thPlot 0])
    hold(a,'all');
    plot(a,[0,0],[-0.5,0.5]*PRF,'w','linewidth',2);
    hold(a,'off');
    xlabel(a,'Doppler frequency [Hz]')
    ylabel(a,'Depth [mm]')
    frameVideo=getframe(h);
    writeVideo(writerObj,frameVideo);
end
close(writerObj);