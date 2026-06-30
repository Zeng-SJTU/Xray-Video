%%%%%%%%%%%%%%%%%%%%%%%%%%% X-ray Video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
path(path,genpath(pwd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Para %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 32; q = 32; key = 1024; threshold = 0.3;
d = 200;                                     % frame rate
rate = 1024;                                 % number of sampling patterns

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Static %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename = '650V1MHz50MSa35kV800uA200Hz0(static)S';  % letter S/J/T/U
% df = 128;                                            % velocity division
% n = 200;                                             % number of stacked cycles
% startfr = 1;                                         % start frame
% endfr = 1;                                           % end frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%velocity division: 128~0.12cm/s; 64~0.24cm/s; 32~0.48cm/s; 16~0.96cm/s; 8~1.92cm/s; 4~3.84cm/s; 2~7.68cm/s;
filename = '650V1MHz50MSa35kV800uA200Hz500(128)';    % 650V1MHz50MSa35kV800uA200Hz500(128/64/32/16/8/4/2)
df = 128;                                            % 128/64/32/16/8/4/2
n = 400;                                             % number of stacked cycles
startfr = 1;                                         % start frame
endfr = 8000;                                        % end frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("CNR\S.mat");                                    % Ground truth of S/J/T/U

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = load([filename,'.mat']);
Syn = str.src1.Data;                                  % Syn signal
Syn = diff(Syn);
Syn(Syn>-0.5) = 0;
Syn = -Syn;
y = str.src2.Data;                                    % bucket signal
[peaksT, locsT] = findpeaks(Syn, 'MinPeakHeight', -3, 'MinPeakDistance', 4000);
T = unique(locsT);

%%%%%%%%%%%%%%%%%%%% Motion Compensatin %%%%%%%%%%%%%%%%%%%%%%%%%%%
CircT = round((T(end)-T(1))/(length(T)-1));
f = 200*2*df;

% %%%%%%%%%%%%%%%%%%%%%%%%% Static %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datalength = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Dynamic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datalength = CircT/f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image reconstruction %%%%%%%%%%%%%%%%%
L = endfr-startfr+1;                  % length of consecutive frames
Video = zeros(32,32,L);
for k = startfr:endfr
    disp(k);
    Ibni = 0;
    St0 = [];
    St1 = [];
    for c = 1:n
        St0(end+1) = T(k+(c-1))+round(datalength*(c-1));
        St1(end+1) = T(k+c)    +round(datalength*c);
    end
    stvec = sort(St0,'ascend');
    spvec = sort(St1,'ascend');
    %%%%%%%%%%%%%%%%%%%%%%%%%% Bucket signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ch2 = 0; vq2 = 0;
    for i0 = 1:length(stvec)
        disp(i0);
        if stvec(i0) > 0
            Ch2 = y(stvec(i0)+1:spvec(i0));
            %%%%%%%%%%%%%%% Interpolation operations %%%%%%%%%%%%%%%%%%%%%%
            DataLens = spvec(i0) - stvec(i0);
            gain = floor(DataLens/key);
            x = linspace(0, 1,  gain*key);
            xq = linspace(0, 1, length(Ch2));
            Ib  = interp1(xq,Ch2,x,'spline','extrap');
            vq2 = Ib + vq2;
            Ch2 = 0;
            spvec(i0) - stvec(i0);
        else
            continue;
        end
    end
    %%%%%%%%%%%%%%%%%%% Evenly distribute %%%%%%%%%%%%%%%%%%%%%%%%%
    ReS = reshape(vq2,[gain key]);
    Reorder = sort(ReS,1);
    ReIb = Reorder(1:end,:);
    vq = fliplr(sum(ReIb,1));
    star = 18;
    Ibn = [vq(star:end) vq(1:star-1) ];
    Ibni = Ibni + Ibn;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = imread('mask10_cu_32x32.bmp');
    size(A);
    Aa = [A;A(1:31,:)];
    for i = 1:key
        S  = double(reshape(Aa(i:i+31,:),[p q]));
        Am(i,:) = double(S(:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TVAL3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear opts
    opts.mu = 2^9;
    opts.beta = 2^6;
    opts.tol = 1E-5;
    opts.maxit = 3000;
    opts.TVnorm = 1;
    opts.nonneg = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ibniv = [Ibni(star:end) Ibni(1:star-1)];
    t1 = cputime;
    [U, out] = TVAL3(Am(1:rate,:),-Ibni(1:rate)', p, q, opts);
    t1 = cputime - t1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BGU2 = CNRbg;
    U_G1 = U(BGU2 == 1);
    U_G0 = U(BGU2 == 0);
    U_CNR = (mean(U_G1(:)) - mean(U_G0(:)))./(sqrt((((std(U_G1(:))).^2) + ((std(U_G0(:))).^2))));
    fprintf('CSGI_CNR：%.2f\n', U_CNR);
    %%%%%%%%%%%%%%%%%%%%%%%%%%% XGI image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U1 = U;
    U1(U1<0) = 0;
    U1 = mat2gray(U1);
    Video(:,:,k-startfr+1) = U1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output.mat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirsave = 'Videos';
outputdata = fullfile(dirsave, [filename, '.mat']);
save(outputdata, 'Video');

%%%%%%%%%%%%%%%%%%%%%%%%%%% Output.video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputVideo = VideoWriter(fullfile(dirsave, [filename, '.avi']), 'Uncompressed AVI');
outputVideo.FrameRate = d;
open(outputVideo);
for fr1 = 1:L
    currentFrame = Video(:,:,fr1);
    currentFrame(currentFrame<threshold) = 0;
    currentFrame = rot90(currentFrame,-1);
    currentFrame = im2uint8(currentFrame);
    writeVideo(outputVideo, currentFrame);
end
close(outputVideo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for fr = 1:L
    U2 = Video(:,:,fr);
    U2(U2<threshold) = 0;
    clf;
    U2 = rot90(U2,-1);
    U2 = flipud(U2);
    U2 = im2uint8(U2);
    U2 = imresize(U2,4,'bilinear');
    rectangle('Position', [1, 1, 138, 148], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    y1 = 3; y2 = 4;
    x_positions = 5:8:134;
    for x = x_positions
        rectangle('Position', [x, y1, 4, y2],'FaceColor', 'k','EdgeColor', 'k','LineWidth', 1);
    end
    y3 = 143; y4 = 4;
    x_positions = 5:8:134;
    for x = x_positions
        rectangle('Position', [x, y3, 4, y4],'FaceColor', 'k','EdgeColor', 'k','LineWidth', 1);
    end
    axis off;
    hold on;
    [rows, cols] = size(U2);
    image([6, 6+127], [11, 11+127], U2, 'CDataMapping', 'scaled');
    colormap('gray');
    hold on;
    NumPixel = (size(U2,2)/128)*16;
    SizePixel = NumPixel*56.25;
    line([5+33,5+96],[10+64,10+64],'Color','red','LineWidth',4,'LineStyle','--');
    line([5+1,5+128],[10+64,10+64],'Color','red','LineWidth',4,'LineStyle','--');
    line([5+64,5+64],[10+1,10+128],'Color','green','LineWidth',4,'LineStyle','--');
    line([116,NumPixel+116],[12,12],'Color','yellow','LineWidth',5);
    text(101+5+NumPixel/2,13,sprintf('%dµm', SizePixel), ...
        'Color','yellow','FontSize',30,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold off;
    t = 1/d;
    pause(t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%