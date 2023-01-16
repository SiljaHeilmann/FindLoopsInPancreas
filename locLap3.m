function [Iout] = locLap3(Iin,sigma,alpha,beta)


%sigma = 0.9;%Amplitude of edges, specified as a non-negative scalar. sigma should be in the range [0, 1] for integer images and for single images defined over the range [0, 1].
%alpha = 0.5;%1.8;%0.05;% Smoothing of details, specified as a positive scalar. Typical values of alpha are in the range [0.01, 10].
% alpha<1 Iinncreases the details of the input image, effectively enhancing the local contrast of the image without affecting edges or introducing halos.
% alpha>1 Smooths details in the input image while preserving crisp edges
%beta =1;% 0.1;%Dynamic range, specified as a non-negative scalar. Typical values of beta are in the range [0, 5]. beta affects the dynamic range of A.
% beta<1 Reduces the amplitude of edges in the image, effectively compressing the dynamic range without affecting details.


locLapZ = double(zeros(size(Iin)));
locLapX = double(zeros(size(Iin)));
locLapY = double(zeros(size(Iin)));


for zz=1:size(locLapZ,3)

locLapZ(:,:,zz) = locallapfilt(uint16(Iin(:,:,zz)),sigma,alpha,beta);

end

for xx=1:size(locLapX,1)

locLapX(xx,:,:) = locallapfilt(uint16(squeeze(Iin(xx,:,:))),sigma,alpha,beta);

end

for yy=1:size(locLapY,2)

locLapY(:,yy,:) = locallapfilt(uint16(squeeze(Iin(:,yy,:))),sigma,alpha,beta);

end

locLap = locLapX;
locLap(locLapY>locLap)=locLapY(locLapY>locLap);
locLap(locLapZ>locLap)=locLapZ(locLapZ>locLap);

Iout = locLap;

end

