function [p_val_matrix,h_matrix,TT_vec,X_tick_positions,Y_median_positions,Y_ci] = boxplotcurve(x_vec,y_vec,num_bins,log_lin_str,color)
% boxplot compact style plot with data divided into num_bins groups and mean x-value for group on x-axis 

if strcmp(log_lin_str,'lin')
    TT_vec = linspace(min(x_vec),max(x_vec),num_bins+1);
    TT_vec(end)=[];
elseif strcmp(log_lin_str,'log')
    TT_vec = logspace(log10(min(x_vec)),log10(max(x_vec)),num_bins+1);
    TT_vec(end)=[];
elseif strcmp(log_lin_str,'Nbin')% same number of points in each bin
    x_vec_sort = sort(x_vec);
    numInOneBin = num_bins;
    count = num_bins;
    TT_vec(1) = x_vec_sort(1);
    cc = 1;
    while count<=length(x_vec)
        cc=cc+1;
        TT_vec(cc)=x_vec_sort(count);% make threshold vec such that each bin has same number of points
        count=count+num_bins;
    end
    num_bins = length(TT_vec)-1;
end


bin_cat_vec = zeros(size(x_vec));
% go through x_vec and decide which number bin each point belongs in
for kk=1:length(x_vec)
    for tt=1:length(TT_vec)-1
      if x_vec(kk)>=TT_vec(tt)%x_vec(kk)>=floor(TT_vec(tt))
         bin_cat_vec(kk)=tt;
      end
    end
end


X_tick_positions = zeros(1,num_bins);
Y_median_positions = zeros(1,num_bins);
Y_cis = zeros(1,num_bins);

for ii=1:num_bins
   X_tick_positions(ii) = median(x_vec(bin_cat_vec==ii));
   X_tick_labels{ii} = ['n=' num2str(length(x_vec(bin_cat_vec==ii)))];
   Y_median_positions(ii) = median(y_vec(bin_cat_vec==ii));
   Y_ci(ii) =  (iqr(y_vec(bin_cat_vec==ii)).*1.57)./sqrt(length(y_vec(bin_cat_vec==ii)));

end



X_tick_positions
Y_median_positions
Y_ci
%unique(bin_cat_vec)
%sum(bin_cat_vec==0)
%TT_vec
%x_vec(bin_cat_vec==0)


boxplot(y_vec,bin_cat_vec,'Positions', X_tick_positions, 'labels',X_tick_labels,'Notch','on','PlotStyle','Compact','Color',color,'symbol','.','jitter',0.1); hold on
if strcmp(log_lin_str,'log')
    set(gca,'xscale','log')
    set(gca,'XTick',[0.1 1 10 100 1000 10000])
    set(gca,'XTickLabel',[0.1 1 10 100 1000 10000])
elseif  strcmp(log_lin_str,'lin')
    set(gca,'XTick',[0 1000 2000 3000 4000])
    set(gca,'XTickLabel',[0 1000 2000 3000 4000])
end

% plot number of points in each bin
text(X_tick_positions,ones(size(X_tick_positions)).*max(y_vec).*(2/3),X_tick_labels)
axis([10.*(floor(min(x_vec-10)./10)) max(x_vec) min(y_vec) max(y_vec)])


p_val_matrix=NaN(num_bins,num_bins);
h_matrix=NaN(num_bins,num_bins);

for ii=1:num_bins
    for jj=ii:num_bins
        if jj~=ii
             [p_val_matrix(ii,jj),h_matrix(ii,jj)] = ranksum(y_vec(bin_cat_vec==ii),y_vec(bin_cat_vec==jj));
        end
    end
end

num_significant = sum(h_matrix);
count=0;
for ii=1:size(h_matrix,1)
    for jj=1:size(h_matrix,2)
        if h_matrix(ii,jj)==1
            count=count+.05;       
            text(mean([X_tick_positions(ii) X_tick_positions(jj)]),min(y_vec)+count,'*')
            plot([X_tick_positions(ii) X_tick_positions(jj)],[min(y_vec)+count min(y_vec)+count],'k-')
        end
    end
end


end

