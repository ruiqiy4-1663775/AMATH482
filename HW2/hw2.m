%% Part 1
clear; close all; clc
load handel  % Fs = the sampling frequency in Hz
             % y = the audio signal amplitude as a single column vector
v = y';
plot((1:length(v))/Fs, v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% p8 = audioplayer(v,Fs);
% playblocking(p8);   

%% different window widths
L = length(v)/Fs;
n = length(v);
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k=(1/L)*[0:(n-1)/2 -(n-1)/2:-1];
ks=fftshift(k);

a = [1 100 1000];
tslide = 0:0.1:L;
vgt_spec = zeros(length(tslide), n);
vgt_spec1 = zeros(length(tslide), n);  
vglt_spec2 = zeros(length(tslide), n);  
for j = 1:length(tslide)
    g = exp(-a(2)*(t-tslide(j)).^2); 
    g_narroww = exp(-a(3)*(t-tslide(j)).^2);
    g_widew = exp(-a(1)*(t-tslide(j)).^2);
    
    vg = g.* v; 
    vgt = fft(vg); 
    vg1 = g_narroww .* v;
    vgt1 = fft(vg1);
    vg2 = g_widew .* v;
    vgt2 = fft(vg2);
    
    vgt_spec(j,:) = fftshift(abs(vgt));
    vgt_spec1(j,:) = fftshift(abs(vgt1));
    vglt_spec2(j,:) = fftshift(abs(vgt2));
end
figure(1)
subplot(1, 3, 1);
pcolor(tslide,ks,vgt_spec.'), 
shading interp
title('a=100')  
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)

subplot(1, 3, 2);
pcolor(tslide,ks,vgt_spec1.'), 
shading interp 
title('Narrow Window: a=1000')
xlabel('Time in sec)'), ylabel('Frequency in Hz')
colormap(hot)

subplot(1, 3, 3);
pcolor(tslide,ks,vglt_spec2.'), 
shading interp 
title('Wide Window: a=1')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)

%% Oversampling and Undersampling
tslide1 = 0:0.01:L;  % increment  = 0.01
vgtspec1 = zeros(length(tslide1), n);
for j = 1:length(tslide1)
    g = exp(-100*(t-tslide1(j)).^2); 
    vg = g.* v;
    vgt = fft(vg);
    vgtspec1(j,:) = fftshift(abs(vgt));
end
tslide2 = 0:1:L;  %  increment = 1
vgtspec2 = zeros(length(tslide2), n);
for j = 1:length(tslide2)
    g = exp(-100*(t-tslide2(j)).^2); 
    vg = g.* v;
    vgt = fft(vg);
    vgtspec2(j,:) = fftshift(abs(vgt));
end
figure(2)
subplot(3,1,1)
pcolor(tslide,ks,vgt_spec.'), 
shading interp 
title('Normal sampling')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)

subplot(3,1,2)
pcolor(tslide1,ks,vgtspec1.'), 
shading interp 
title('Oversampling')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)

subplot(3,1,3)
pcolor(tslide2,ks,vgtspec2.'), 
shading interp 
title('Undersampling')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)


%% Part I Using Different Gabor windows

vmt_spec = zeros(length(tslide), n); 
vst_spec = zeros(length(tslide), n);
for j = 1:length(tslide)
    %  Mexican hat wavelet
    sigma = 0.05;
    mex_hat = (2/(sqrt(3*sigma)*(pi^0.25))).*(1-((t-tslide(j))/sigma).^2)...
        .* exp(-((t-tslide(j)).^2)/(2*sigma^2));

    % Step-function (Shannon) window
    width = 0.05;
    shn = abs(t - tslide(j)) <= width /2;
    
    vmh = mex_hat .* v;
    vmht = fft(vmh);
    
    vsh = shn .* v;
    vsht = fft(vsh);

    vmt_spec(j,:) = fftshift(abs(vmht));
    vst_spec(j,:) = fftshift(abs(vsht));
end

figure(3)
subplot(3,1,1)
pcolor(tslide,ks,vgt_spec.'), 
shading interp 
title('Gaussian function')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)

subplot(3,1,2)
pcolor(tslide,ks,vmt_spec.'), 
shading interp 
title('Mexican hat function')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)

subplot(3,1,3)
pcolor(tslide,ks,vst_spec.'), 
shading interp 
title('Step-function')
xlabel('Time in sec'), ylabel('Frequency in Hz')
colormap(hot)


%% Part II - Starter Code
clear; close all; clc

[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');
% p8 = audioplayer(y,Fs); playblocking(p8);


%% Part II - 1
% piano
v = y';
n = length(v);
t2 = linspace(0,tr_piano,n+1); 
t_p = t2(1:n); 
k_p = (2*pi/tr_piano)*[0:n/2-1 -n/2:-1]; 
ks_p = fftshift(k_p);

tslide_p = 0:0.2:tr_piano;
spectrogram_p = zeros(length(tslide_p), n);
notes_piano = zeros(1, length(tslide_p));
for j = 1:length(tslide_p)
    f = exp(-100*(t_p-tslide_p(j)).^2);
    window = f .* v;
    windowfft = fft(window);
    [M, I] = max(windowfft); %Find max frequency
    notes_piano(1,j) = abs(k_p(I))/(2*pi);
    spectrogram_p(j,:) = fftshift(abs(windowfft));
end

% recorder
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (recorder)');
% p8 = audioplayer(y,Fs); playblocking(p8);

v = y';
n = length(v);
t2 = linspace(0,tr_rec,n+1); 
t_rec = t2(1:n); 
k_rec = (2*pi/tr_rec)*[0:n/2-1 -n/2:-1]; 
ks_rec = fftshift(k_rec);

tslide_rec = 0:0.2:tr_rec;
spectrogram_rec = zeros(length(tslide_rec), n);
notes_rec = zeros(1, length(tslide_rec));
for j = 1:length(tslide_rec)
    f = exp(-100*(t_rec-tslide_rec(j)).^2);
    
    window = f .* v;
    windowfft = fft(window);
    [M, I] = max(windowfft);  %Find max frequency
    notes_rec(1,j) = abs(k_rec(I))/(2*pi);
    spectrogram_rec(j,:) = fftshift(abs(windowfft));
end


figure(4)
subplot(2,1,1)
pcolor (tslide_p,(ks_p/(2*pi)), spectrogram_p .'), 
shading interp
title ("Piano")
xlabel('Time in sec'), ylabel('Frequency in Hz')
ylim ([0 400])
colormap(hot)

subplot(2,1,2)
pcolor (tslide_rec,(ks_rec/(2*pi)), spectrogram_rec .'), 
shading interp
title ("Recorder")
xlabel('Time in sec'), ylabel('Frequency in Hz')
ylim ([0 1400])
colormap(hot)

%%
% Reproduce the music score/note on the piano and recorder
figure(5)
subplot(2, 1, 1)
plot(tslide_p, notes_piano,'o','MarkerFaceColor', 'b');
yticks([246.9417,261.6256,293.6648,329.6276,349.2282]); 
yticklabels({'B4','C4','D4','E4','F4'});
ylim ([246 350])
title("Score for Piano Music");
xlabel("Time (s)"); ylabel("Notes corresponding to frequency (Hz)");

subplot(2, 1, 2)
plot(tslide_rec, notes_rec,'o','MarkerFaceColor', 'b');
yticks([783.99, 880, 987.77, 1046.5, 1174.7, ]); 
yticklabels({'G5','A6','B6','C6','D6'});
ylim ([700 1200])
title("Score for Recording");
xlabel("Time (s)"); ylabel("Notes corresponding to frequency (Hz)");




