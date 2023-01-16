function [im_eq,h_eq,h_im,h_ref,T_map] = nanhisteq(im,im_ref,nbins,show_yes_no)
% returns an version of im which now has a histogram that aproximately
% matches the histogram om im_ref
% input images should be type double and they may contain NaN's

 if strcmp('double',class(im))==0 || strcmp('double',class(im_ref))==0 
    display('Silja: Error im and im_ref need to be type double')
 end
 if max(im_ref(:))>1 || max(im(:))>1 || max(im_ref(:))<0 || max(im(:))<0
    display('Silja: Error im and im_ref need to be in the range [0 1]')
    max_im = max(im(:))
    max_im_ref = max(im_ref(:))
 end
 
 edges = linspace(0,1,nbins);
 
 h_ref = histcounts(im_ref(isfinite(im_ref)),edges,'Normalization','Probability');
 
 [im_eq_no_nan,T] = histeq(im(isfinite(im)),h_ref);% T is a 256 long grayscale tranformation map
 
 h_eq = histcounts(im_eq_no_nan,edges,'Normalization','Probability');
 
 im_eq = zeros(size(im));
 
 index = linspace(0,1,256);
 
 for tt=2:length(T)
     temp = im;
     temp(im>=index(tt))=0;
     temp(im<index(tt-1))=0;
     
     im_eq(temp>0)=T(tt);     

%      temp = im;
%      temp(im>=T(tt))=0;
%      temp(im<T(tt-1))=0;
%      
%      im_eq(temp>0)=index(tt); 
 end
 
 im_eq(isfinite(im)==0) = NaN;
 
 h_eq2 = histcounts(im_eq(isfinite(im_eq)),edges,'Normalization','Probability');
 
 h_im = histcounts(im(isfinite(im)),edges,'Normalization','Probability');

 
 switch nargin
    case 3
        disp('Silja: put in a random 4th argument to get grafical output of image histograms.')
    case 4

        
        figure
        subplot(231); hold on
        imshow(im, [0 0.4]); colorbar; title('im')
        subplot(232); hold on
        imshow(im_ref, [0 0.4]); colorbar; title('im-ref')
        subplot(233); hold on
        imshow(im_eq, [0 0.4]); colorbar; title('im-eq')
        subplot(234); hold on
        plot(edges(1:end-1),h_im,'.'); set(gca,'yscale','log','xscale','log'); legend('im')%,'xscale','log'
        subplot(235); hold on
        plot(edges(1:end-1),h_im,'.')
        plot(edges(1:end-1),h_ref,'.');set(gca,'yscale','log','xscale','log'); legend('im','im-ref')
        subplot(236); hold on
        plot(edges(1:end-1),h_eq,'.')
        plot(edges(1:end-1),h_ref,'.');set(gca,'yscale','log','xscale','log'); legend('im-eq','im-ref')
        plot(edges(1:end-1),h_eq2,'.')
        
%                 figure
%         subplot(231); hold on
%         imshow(im, [0 0.4]); colorbar; title('im')
%         subplot(232); hold on
%         imshow(im_ref, [0 0.4]); colorbar; title('im-ref')
%         subplot(233); hold on
%         imshow(im_eq, [0 0.4]); colorbar; title('im-eq')
%         subplot(234); hold on
%         plot(edges(1:end-1),cumsum(h_im),'.'); set(gca,'yscale','log','xscale','log'); legend('im')%,'xscale','log'
%         subplot(235); hold on
%         plot(edges(1:end-1),cumsum(h_im),'.')
%         plot(edges(1:end-1),cumsum(h_ref),'.');set(gca,'yscale','log','xscale','log'); legend('im','im-ref')
%         subplot(236); hold on
%         plot(edges(1:end-1),cumsum(h_eq),'.')
%         plot(edges(1:end-1),cumsum(h_ref),'.');set(gca,'yscale','log','xscale','log'); legend('im-eq','im-ref')
%         plot(edges(1:end-1),cumsum(h_eq2),'.')
        
        figure
        plot(index,T,'.'); hold on
        xlabel('index')
        ylabel('T map')
        plot(index,index,'-')
        
        axis([0 max(T)-1/56 0 max(T)]-1/56); axis square
        
    otherwise
        disp('Silja: wrong number of arguments in')
end
    
T_map =T;
 
end

