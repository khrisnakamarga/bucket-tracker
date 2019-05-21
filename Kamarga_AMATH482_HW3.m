% AMATH 482 - HW3
% Khrisna Kamarga

% Setup - Load the Video and Prepare Variables
clear all; close all; clc;
load cam1_4.mat;
video = (vidFrames1_4);
[m, n, rgb, t] = size(video);

%% visualize the motion for each dataset
for i = 140:t
%     bwVideo = rgb2gray(video(:,:,:,i));
%     pcolor((bwVideo)), shading interp, colormap(gray);
    imshow(video(:,:,:,i)); drawnow;
%     drawnow
end

%% crop location for 1_1 & 1_2 & 1_3 & 1_4
isolate = zeros(m,n);
jStart = 250;
iStart = 150;
jEnd = 450;
iEnd = 450;
%% 2_1 & 2_2 & 2_3 & 2_4
isolate = zeros(m,n);
jStart = 150;
iStart = 100;
jEnd = 400;
iEnd = 400;
%% 3_1 & 3_2 & 3_3
isolate = zeros(m,n);
jStart = 200;
iStart = 200;
jEnd = 500;
iEnd = 350;

%% 3_4
isolate = zeros(m,n);
jStart = 300;
iStart = 150;
jEnd = 500;
iEnd = 300;

%% create the crop filter
for i = iStart:iEnd
    for j = jStart:jEnd
        isolate(i,j) = 1;
    end
end
% crop filter preview
clc;
% for i = 1:50
%     bwVideo = rgb2gray(video(:,:,:,i));
%     pcolor(flipud(double(bwVideo).*isolate)), shading interp, colormap(gray);
%     drawnow
% end

% find characteristic frequency
% set up the fourier coefficients
xaxis2 = linspace(0,m,m+1); x = xaxis2(1:m);
yaxis2 = linspace(0,n,n+1); y = yaxis2(1:n);
kx = (2*pi/(2*m))*[0:(m/2-1) -m/2:-1]; ksx = fftshift(kx);
ky = (2*pi/(2*n))*[0:(n/2-1) -n/2:-1]; ksy = fftshift(ky);
% set up the 3D coordinate points
[X,Y] = meshgrid(x,y); % spatial coordinates
[Kx,Ky] = meshgrid(ksx,ksy); % wave numbers
% transposing the coordinates since the pictures cause weird flips
X = X';
Y = Y';
Kx = Kx';
Ky = Ky';

UtnAve = zeros(m,n); % kernel for the averaged frequency domain signal
background = zeros(m,n); % background kernel
for i = 1:t
    bwVideo = double(rgb2gray(video(:,:,:,i)));
    Un = bwVideo; % gets the 2D coordinate representation of the sample
    background = background + double(Un); 
    Utn = fftn(Un); % fourier transform of the data
    UtnAve = UtnAve + Utn; % cummulative sum of the frequency domain signal
end

%average background
background = background / t;

% %% plot the spectrogram
% % plot the resulting normalized averaged data in the frequency domain
% pcolor(abs(fftshift(UtnAve))/max(abs(UtnAve), [], 'all')), shading interp, colormap(gray);
% xlabel("Kx"); ylabel("Ky");
% title("Averaged Data in the Frequency Domain");

% find the indices of the max magnitude in the frequency domain
[ind1 ind2] = ind2sub([m,n], find(fftshift(UtnAve) == max(fftshift(UtnAve), [], 'all')));
% look up the frequency domain coordinate of the strongest signal
Kc = [Kx(ind1, ind2), Ky(ind1, ind2)];

% applying the filter and recovering x and y
close all; clc;
% 2D gaussian filter
tau = 100; % bandwith of the filter (good: 0.2)
[kux, kuy] = meshgrid(kx,ky); % unshifted wave numbers
kux = kux';
kuy = kuy';
filter = exp(-tau*((kux - Kc(1)).^2+(kuy - Kc(2)).^2));

bucket = zeros(t, 2); % kernel for the coordinates of the bucket
for i = 1:t
    bwVideo = double(rgb2gray(video(:,:,:,i)));
    Un = (bwVideo - background).*isolate; % gets the 2D coordinate representation of the sample
    Utn = fftn(Un); %Utn = fftshift(Utn);
    UtnFilter = Utn.*filter; % filtered frequency domain signal
    UnFilter = real(ifftn(UtnFilter)); % obtain the spatial filtered data
    
    % draw the resulting spatial filtered data
    pcolor(flipud(abs(UnFilter)/max(abs(UnFilter), [], 'all')))
    shading interp, colormap(gray); 
    grid on

    % find the coordinate of the center of the bucket
    [ind1 ind2] = ind2sub([m,n], find(abs(UnFilter) == max(abs(UnFilter), [], 'all')));
    bucket(i,:) = [X(ind1, ind2), Y(ind1, ind2)];
    
    hold on
    plot(bucket(i,2), -(bucket(i,1)-m), 'm.-', 'MarkerSize', 20);
    hold off
    drawnow
end
% %% plot the trajectory of the bucket
% plot(bucket(:,2), bucket(:,1), 'm.-', 'MarkerSize', 20);
% xlim([0 m]); ylim([0 n]);
% title("Trajectory of the bucket");
% xlabel("X"); ylabel("Y");
% grid on
% %% plot the decoupled x and y data
% subplot(2,1,1)
% plot(1:length(bucket(:,1)), bucket(:,1));
% title("X");
% subplot(2,1,2)
% plot(1:length(bucket(:,2)), bucket(:,2));
% title("Y");

% filter the noisy x and y
% set up the fourier coefficients
xaxis2 = linspace(0,t,t+1); xaxis = xaxis2(1:t);
k = [0:(t/2-1) -t/2:-1]; ks = fftshift(k);
% 3D gaussian filter
tau = 0.001; % bandwith of the filter (good: 0.2)
filter = exp(-tau*(k).^2);
filter = filter';

x = real(ifft(fft(bucket(1:length(k),2)).*filter));
y = real(ifft(fft(bucket(1:length(k),1)).*filter));

subplot(2,1,2)
plot(1:length(x), x)
title("X-axis displacement")
subplot(2,1,1)
plot(1:length(y), y)
title("Y-axis displacement")

clearvars -except x y

%% save variables
x1=x;
y1=y;
save('vid1')
%%
x2=x;
y2=y;
save('vid2')
%%
x3=x;
y3=y;
save('vid3')
%% PCA
clear all; clc; close all;
load('vid1.mat')
load('vid2.mat')
load('vid3.mat')
clearvars -except x1 x2 x3 y1 y2 y3
x1 = x1';
x2 = x2';
y2 = y2';

% maxSize = max([length(x1) length(x2) length(x3) length(y1) length(y2) length(y3)]);
maxSize = 226;
X = zeros(6,maxSize);
X(1,:) = x1(1:226);
X(2,:) = y1(1:226);
X(3,:) = x2(1:226);
X(4,:) = y2(1:226);
X(5,:) = x3(1:226);
X(6,:) = y3(1:226);

%%
clear all; clc;
load '4.mat'
[m,n]=size(X); % compute data size
mn=mean(X,2); % compute mean for each row
X=X-repmat(mn,1,n); % subtract mean
Cx=(1/(n-1))*X*X'; % covariance
[u,s,v]=svd(X/sqrt(n-1)); % perform the SVD
lambda=diag(s).^2; % produce diagonal variances
Y=u'*X; % produce the principal components projection

% figure(1)
% for j = 1:6
%     subplot(3,2,j);
%     plot(1:226,-Y(j,1:226));
%     title("Principal Component " + j);
% end

%%
hold on
figure(2)
plot(1:length(lambda), lambda, 'bo','LineWidth',3)

%%
% save('4')

%% Results
clear all; close all; clc;
load '3.mat'

%%
close all; clc;
figure(1)
for j = 1:6
    subplot(3,2,j);
    plot(1:226,Y(j,1:226));
    title("Principal Component " + j);
    xlabel("index");
end

figure(2)
plot(1:length(lambda), lambda);
title("Energy Plot")
xlabel("Principal Component")
ylabel("Energy")

figure(3)
Cy = cov((u'*X)');
pcolor(flipud(Cy)), shading interp, colormap
title("Energy Matrix (S)")
ylabel("row");
xlabel("column");



