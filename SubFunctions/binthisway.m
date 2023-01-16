function [p_val_matrix,h_matrix] = binthisway(x_vec,y_vec,TT_vec,color,linestyle)
% plot with medians of binned data with data divided into length(TT_vec) groups and mean x-value for group on x-axis 
% and confidence intervals as errorbars


bin_cat_vec = zeros(size(x_vec));
% fordel data i bins
for kk=1:length(x_vec)
    for tt=1:length(TT_vec)
      if x_vec(kk)>=TT_vec(tt)%floor(TT_vec(tt))
         bin_cat_vec(kk)=tt;
      end
    end
end

X_tick_positions = zeros(1,length(TT_vec));

% find center for each data bin
for ii=1:length(TT_vec)
   X_tick_positions(ii) = mean(x_vec(bin_cat_vec==ii));
   X_tick_labels{ii} = ['n=' num2str(length(x_vec(bin_cat_vec==ii)))];
end

y_vec_plot=zeros(1,length(X_tick_positions));
error_vec=zeros(1,length(X_tick_positions));

for ii=1:length(X_tick_positions)
    y_vec_plot(ii) = median(y_vec(bin_cat_vec==ii)); 
    prc = prctile(y_vec(bin_cat_vec==ii),[25 50 75]);
    error_vec(ii) = 1.57*(prc(3)-prc(1))/sqrt(length(y_vec(bin_cat_vec==ii)));% 0.5 length og error bar
    % 1.57*(Q75-Q25)/sqrt(length(data)) How to calculate CI's
end

% X_tick_positions

% y_vec_plot
% error_vec
% vec=unique(bin_cat_vec);
% for ii=1:length(unique(bin_cat_vec))
%     vec(ii)
%     sum(bin_cat_vec==vec(ii))
% end

errorbar(X_tick_positions,y_vec_plot,error_vec,'color', color,'Linewidth',2,'Linestyle',linestyle); hold on
%plot(X_tick_positions,y_vec_plot,'.','color', color,'MarkerSize',10); 

% plot number of points in each bin
%text(X_tick_positions,ones(size(X_tick_positions)).*max(y_vec).*(1/4),X_tick_labels)
text(X_tick_positions-50,y_vec_plot,X_tick_labels)

axis([10.*(floor(min(x_vec-10)./10)) max(x_vec) min(y_vec) max(y_vec)])


p_val_matrix = NaN(size(TT_vec,2),size(TT_vec,2));
h_matrix = NaN(size(TT_vec,2),size(TT_vec,2));


for ii=1:size(TT_vec,2)
    for jj=ii:size(TT_vec,2)
        if jj~=ii
            if isempty(find(isfinite(y_vec(bin_cat_vec==ii)),1))==0 && isempty(find(isfinite(y_vec(bin_cat_vec==jj)),1))==0 
                [p_val_matrix(ii,jj),h_matrix(ii,jj)] = ranksum(y_vec(bin_cat_vec==ii),y_vec(bin_cat_vec==jj));
            end
        end
    end
end

f=0;%2.95;
count=0;
delta=0;%.005;
for ii=1:size(h_matrix,1)
    for jj=1:size(h_matrix,2)
        if p_val_matrix(ii,jj)<=0.0001
            count=count-delta;       
            text(mean([X_tick_positions(ii) X_tick_positions(jj)]),min(y_vec)*f+count,'****')
            plot([X_tick_positions(ii) X_tick_positions(jj)],[min(y_vec)*f+count min(y_vec)*f+count],'k-')
        elseif p_val_matrix(ii,jj)<=0.001
            count=count-delta;       
            text(mean([X_tick_positions(ii) X_tick_positions(jj)]),min(y_vec)*f+count,'***')
            plot([X_tick_positions(ii) X_tick_positions(jj)],[min(y_vec)*f+count min(y_vec)*f+count],'k-')
        elseif p_val_matrix(ii,jj)<=0.01
            count=count-delta;       
            text(mean([X_tick_positions(ii) X_tick_positions(jj)]),min(y_vec)*f+count,'**')
            plot([X_tick_positions(ii) X_tick_positions(jj)],[min(y_vec)*f+count min(y_vec)*f+count],'k-')
        elseif p_val_matrix(ii,jj)<=0.05
            count=count-delta;       
            text(mean([X_tick_positions(ii) X_tick_positions(jj)]),min(y_vec)*f+count,'*')
            plot([X_tick_positions(ii) X_tick_positions(jj)],[min(y_vec)*f+count min(y_vec)*f+count],'k-')
        elseif p_val_matrix(ii,jj)>0.05
            count=count-delta;       
            text(mean([X_tick_positions(ii) X_tick_positions(jj)]),min(y_vec)*f+count,'ns')
            plot([X_tick_positions(ii) X_tick_positions(jj)],[min(y_vec)*f+count min(y_vec)*f+count],'k-')       
        end
    end
end

end




