function [peakL,Hlist,Hrel,i_found,acc] = eval_fit(peaks,peaks_fit,H,H_fit,delta,verbose,label,fnam,n_assign)
%Evaluates fit to synthetic data
%   Detailed explanation goes here

addTitle = false;
mkFile = false;
compF = false;
doH = true;
minDiff = 0.005;

if nargin > 4 && ~isempty(delta)
    minDiff = delta;
end
if nargin < 6
    verbose = 1;
elseif nargin > 6 && ~isempty(label)
    addTitle = true;
elseif nargin > 7 && ~isempty(fnam)
    mkFile = true;
end

% if isempty(H_fit)
%     compF = true;
%     doH = false;
%     minDiff = 0.0001;
% end
% 
% 
% if ~(size(peaks,1) == size(peaks_fit,1))
%     % Make peaklists comparable length by removing integer masses not in
%     % both lists. Assumes one of the lists contains all integer masses in 
%     % other list.
%     if (size(peaks,1) < size(peaks_fit,1))
%         ind = any(round(peaks_fit(:,1)) == round(peaks(:,1))',2);
%         peaks_fit = peaks_fit(ind,:);
%     else
%         ind = any(round(peaks(:,1)) == round(peaks_fit(:,1))',2);
%         peaks = peaks(ind,:);
%     end
% end

numbers = sum(~isnan(peaks),2);
numbers_fit = sum(~isnan(peaks_fit),2);
n = length(numbers);
rate = sum(numbers==numbers_fit)/n;
rate_err = sqrt((rate*(1-rate))/n);

peakL = nan(sum(~isnan(peaks),'all')+sum(~isnan(peaks_fit),'all'),2);
Hlist = peakL;
Hrel = Hlist;
diffvec = Hrel(:,1);
j = 1;

if length(minDiff) == 1
    minDiff = ones(n,1).*minDiff;
end


for i = 1:n
    [peaks1,peaks2,ind1,ind2] = custom_sort(peaks(i,:),peaks_fit(i,:),minDiff(i));
    peakn = length(peaks1);
    peakL(j:j+peakn-1,1) = peaks1;
    peakL(j:j+peakn-1,2) = peaks2;
    diffvec(j:j+peakn-1) = minDiff(i);

    Hlist(j+ind1-1,1) = sum(H(1:sum(~isnan(peaks(i,:))),:,i),2);
    
    if doH
        Hlist(j+ind2-1,2) = sum(H_fit(1:sum(~isnan(peaks_fit(i,:))),:,i),2);
    end
    Hrel(j+ind1-1,1) = Hlist(j+ind1-1,1)/sum(Hlist(j+ind1-1,1));
    Hrel(j+ind1-1,2) = Hlist(j+ind1-1,2)/sum(Hlist(j+ind1-1,2));
    
    j = j+peakn;
end

% peakL(:,1): peaks generated
% peakL(:,2): peaks fitted

peakL(j:end,:) = [];
Hlist(j:end,:) = [];
diffvec(j:end,:) = [];

diff = abs(peakL(:,2)-peakL(:,1));
n_wrong = sum(diff > diffvec,'omitnan'); % Mistakenly coupled peaks
diff(diff > diffvec) = nan;

n = sum(~isnan(peaks),'all');
miss_rate = (sum(isnan(peakL(:,2)))+n_wrong)/n;
miss_err = sqrt(miss_rate*(1-miss_rate)/n);

n = sum(~isnan(peaks_fit),'all');
true_rate = 1-(sum(isnan(peakL(:,1)))+n_wrong)/n;
true_err = sqrt(miss_rate*(1-true_rate)/n);

if compF
    n = n_assign;
    assign_rate = 1-(sum(isnan(peakL(:,1)))+n_wrong)/n;
    assign_err = sqrt(miss_rate*(1-true_rate)/n);
end

[mass_rateGen, mass_rateFit] = corrFitFraction(peakL(:,1),peakL(:,2),Hlist(:,1),Hlist(:,2));


%acc:
% 1. Fraction where number of peaks is same as number generated (useless)
% 2. Fraction of generated peaks found
% 3. Fraction of fitted peaks correct
% 4. Fraction of signal identified
% 5. fraction of signal correctly fit (practically identical with 4)
acc = [rate, 1-miss_rate,true_rate,mass_rateGen,mass_rateFit];
% [fitnumber == gennumber, found, correct, signal found, signal correct]

if addTitle
    title_str = "Fit accuracy distribution, "+label;
else
    title_str = "Fit accuracy distribution";
end

if verbose
    fprintf('Fraction of correct number of peaks: %4.2f \x00B1 %4.2f\n',rate,rate_err)
    fprintf('Fraction of generated peaks found: %4.2f \x00B1 %4.2f\n',1-miss_rate,miss_err)
    fprintf('Fraction of fitted peaks correct: %4.2f \x00B1 %4.2f\n',true_rate,true_err)
    fprintf('Fraction of signal found: %4.2f',mass_rateGen)
    fprintf('Fraction of signal correct: %4.2f',mass_rateFit)
    if compF
        fprintf('Fraction of assigned peaks correct: %4.2f \x00B1 %4.2f\n',assign_rate,assign_err)
    end
    
    figure;
    histogram(diff,0:0.001:0.01)
    %title(title_str)
    xlabel("Distance between generated peak and fit \Delta m/z")
    ylabel("Number of peaks")
    xlim([0 0.01])
end

i_found = ~isnan(peakL(:,1)) & ~isnan(peakL(:,2)) & (abs(peakL(:,2)-peakL(:,1))<=0.01);
dist = abs(peakL(:,1)-peakL(:,1)');
dist(dist == 0) = nan;
mind = min(dist,[],2);
mind(mind > 0.5) = nan;

if verbose
    figure;
    scatter(mind(i_found),Hrel(i_found,1),'filled');
    hold on
    scatter(mind(~i_found),Hrel(~i_found,1),'rx','LineWidth',2);
    hold off
    xlabel("\Delta m/z to closest fit")
    ylabel("Relative intensity")
    legend(["Hittad","Missad"])
    
    figure;
    scatter(mind(i_found),Hlist(i_found,1),'filled');
    hold on
    scatter(mind(~i_found),Hlist(~i_found,1),'rx','LineWidth',2);
    hold off
    xlabel("\Delta m/z to closest fit")
    ylabel("Absolute peak size")
    legend(["Found","Missed"])
end

if mkFile
    makeTestTable(peakL,Hlist,fnam);
end

end

function [fracGen, fracFit] = corrFitFraction(peaks,fits,Hgen,Hfit)
    is_corr = ~isnan(peaks) & ~isnan(fits);
    %H_corr = sum(min(sum(Hgen(is_corr,:),2),sum(Hfit(is_corr,:),2)),'all','omitnan');
    H_corr = sum(min(Hgen(is_corr),Hfit(is_corr)),1,'omitnan');
    H_totGen = sum(Hgen,'all','omitnan');
    H_totFit = sum(Hfit,'all','omitnan');
    fracGen = H_corr/H_totGen;
    fracFit = H_corr/H_totFit;
end
