function [BW] = otsu3D(IM,factor)
% finds otsu's thres T and returns binarized image with T*factor as
% threshold

IM(isinf(IM))=0;

if 0>min(IM(:))
    display('A pixel value was below zero! IM = IM-min(IM(:));')
    IM = IM-min(IM(:));
end

IM = IM-min(IM(:)); % make sure image minimum is zero

T = otsuthresh(IM(:)).*max(IM(:));


% 
% figure; hold on
% h=histogram(IM(:));
% xlabel('pixel intensities')
% ylabel('abundance')
% set(gca,'YScale','log')
% plot([T*factor T*factor],[0.000001 max(h.Values)],'-r','LineWidth',2)
% title(['Otsus threshold multiplied by ' num2str(factor)]) 

T*factor

BW = IM>T*factor;
end

