function ACCENSE_cluster(y_out,X,X_all,markers,nfolder,sig_opt)
  
close all 
%arg-check
if nargin < 2
    flag_X = 0;
else
    flag_X = 1;
    if size(y_out,1) ~= size(X,1)
    disp('ERROR : y_out and X should have the same number of datapoints. Exiting. Please check input...');
    return
    end
end

%Assign markers numbers if the variable does not exist
if ~exist('markers','var') && flag_X == 1
    disp('Marker names have not been provided. Using integers 1, 2, ...');
    for i=1:size(X,2)
        markers{i} = num2str(i);
    end
end 
 
% Create folders
mkdir([nfolder '\clusters'])
mkdir([nfolder '\barplots'])

%Find optimal kernel bandwidth
if ~exist('sig_opt','var')
    disp('Finding the optimal kernel bandwidth. This might take a couple of hours (consider parallelizing).. ')
    y_range = min([max(y_out(:,1))-min(y_out(:,1)); max(y_out(:,2)) - min(y_out(:,2))]);
 
    sig_min = y_range/200; sig_max = y_range/20;
    sig_range = linspace(sig_min,sig_max,30);
  
    parfor i=1:length(sig_range) %Change 'for' to 'parfor' if you wish to parallelize
 
        fprintf('Evaluating clusters at bandwidth = %3.1f\n',sig_range(i));
        [SNEmap, ~, ~] = kernel_density(y_out(:,1),y_out(:,2),sig_range(i));
        p = DetectPeaks(flipud(SNEmap));
        Npeaks(i) = size(p,1)/2;
 
    end
 
    %Locate plateau
    for i=1:length(sig_range)
        if i > 1
            if Npeaks(i-1) - Npeaks(i) <= 1
                sig_opt = sig_range(i-1);
                Nsub = Npeaks(i-1);
                break;
            end
        end
    end
 
    figure
    plot(sig_range,Npeaks,'k','LineWidth',3);
    set(gca,'FontSize',14);
    xlabel('$\gamma$','Interpreter','Latex','FontSize',20);
    ylabel('$N_{peaks}$','Interpreter','Latex','FontSize',20);
    saveas(gcf, [nfolder '\Npeaks_bandwidth'],'png'); hold off;
 
    
    if ~exist('sig_opt','var')
        disp('Could not locate a plateau. Consider changing the kernel search parameters on line 50-51');
        reply = input('Would you like to provide a bandwidth (y/n):', 's');
        if reply == 'y'
            sig_opt = input('Please provide a value (Warning : This will directly determine # of subpopulations):');
        else
            disp('Exiting. Consider changing kernel search parameters');
            return
        end
    else
        fprintf('Based on plateau in Npeaks vs. gamma, the optimal kernel bandwidth selected is %3.2f,\n and %i subpopulations are identified\n',sig_opt,Nsub);
        reply = input('Do you want to overrule this and provide your own value of the bandwidth? (y/n):', 's');
        disp(sig_opt)
        if reply == 'y'
            sig_opt = input('Please specify the desired kernel bandwidth (Warning : This will directly determine # of subpopulations):');
        end
    end
end
%Compute density map
[SNEmap, z1_range, z2_range] = kernel_density(y_out(:,1),y_out(:,2),sig_opt);
 
%Find peaks (Returns pixel coordinates)
p = DetectPeaks(flipud(SNEmap));
 
%Convert pixel coordinates to actual coordinates
peaks = [z1_range(p(2:2:end))' z2_range(length(z2_range) - p(1:2:end))'];

% scatter plot
l=length(y_out);
scatter(y_out(:,1),y_out(:,2),ones(l,1)*0.5,'filled'); hold on
alpha(0.4)
L=max(max(abs(y_out)));
xlim([-L L])
ylim([-L L])
xlabel('tSNE_1','FontSize',12);
ylabel('tSNE_2','FontSize',12);
title('Scatter Plot');
axis equal
saveas(gcf,[nfolder '\clusters\Scatter Plot'],'tiff'); hold off;
 
%Show 2-d density maps with subpopulation locations
figure
imagesc(flipud(SNEmap)); hold on;
plot(p(2:2:end),p(1:2:end),'ko','MarkerSize',5,'MarkerFaceColor','k')
colormap(flip(hot));
xlabel('$y_1$','Interpreter','Latex','FontSize',20);
ylabel('$y_2$','Interpreter','Latex','FontSize',20);
title('Subpopulation Centers','FontSize',20);
set(gca, 'XTick', [1 400]);
set(gca, 'XTickLabel', [min(y_out(:,1)) max(y_out(:,1))],'FontSize',12);
set(gca, 'YTick', [1 400]);
set(gca, 'YTickLabel', fliplr([min(y_out(:,2)) max(y_out(:,2))]),'FontSize',12);
colorbar
saveas(gcf,[nfolder '\subpopulations_marked'],'png'); hold off;
 
 
 
%Number subpopulations
figure
imagesc(flipud(SNEmap)); hold on;
colormap(flip(hot))
for i=1:size(p,1)/2
    text(p(2*i,1),p(2*i-1,1),num2str(i),'FontSize',14,'FontWeight','bold','Color','k');
end
xlabel('$y_1$','Interpreter','Latex','FontSize',20);
ylabel('$y_2$','Interpreter','Latex','FontSize',20);
set(gca, 'XTick', [1 400]);
set(gca, 'XTickLabel', [min(y_out(:,1)) max(y_out(:,1))],'FontSize',12);
set(gca, 'YTick', [1 400]);
set(gca, 'YTickLabel', fliplr([min(y_out(:,2)) max(y_out(:,2))]),'FontSize',12);
colorbar
saveas(gcf,[nfolder '\Numbered_subpopulations'],'png'); hold off;
 
 
%If X is not present, return
 
if flag_X == 0
    return
end
 
figure
%Obtain subpopulation samples
for i=1:size(peaks,1)
    
    %Find closest subpopulation and determine appropriate radius to sample
    dists = sqrt(sum((repmat(peaks(i,:),size(peaks,1)-1,1) - [peaks(1:i-1,:); peaks(i+1:size(peaks,1),:)]).^2,2));
    
    min_dist = min(dists);
    
    ind = find(sum((y_out - repmat(peaks(i,:),size(y_out,1),1)).^2,2) < (min_dist/2)^2);
    
    Subpop{i}.X = X(ind,:);
    Subpop{i}.size = size(Subpop{i}.X,1);
    Subpop{i}.X_all = X_all(ind,:);
    Subpop{i}.Ind=ind;
    
    
    scatter(y_out(ind,1),y_out(ind,2),0.5,'filled'); hold on;
    alpha(0.4)
    xlim([-L L])
    ylim([-L L])
    xlabel('tSNE_1','FontSize',12);
    ylabel('tSNE_2','FontSize',12);
    title('Cells employed in subpopulation samples to compute phenotypic signatures');
    saveas(gcf,[nfolder '\clusters\Cluster_' num2str(i)],'png'); hold off;
    hold off;
    
end
save([nfolder '\sub.mat'],'Subpop','sig_opt')
 
 
disp('Computing the Phenoptyic signatures of all subpopulations and storing them in a directory PhenSig\n');
 
 
med_all = median(X);
 
 
for niche = 1:length(Subpop)
    
    cell_data = Subpop{niche}.X;
    bar(median(cell_data),'b');hold on;
    alpha(0.4)
    xlabel('Features','FontSize',15,'FontWeight','b');
    ylabel('Signal','FontSize',15,'FontWeight','b');
    set(gcf,'units','inches');
    set(gcf,'papersize',[12,8])
    set(gcf,'Position',[0 0 10.12 6.56]);
    errorbar([1:size(X,2)],median(cell_data,1),mad(cell_data,[],1)/2,'.k');
    
    %Plot median of full population
    for j=1:length(markers)
        plot([j-0.4 j+0.4],[med_all(j) med_all(j)],'LineWidth',2); hold on;
    end
    
    title(['Subpopulation ' num2str(niche) ', N = ' num2str(size(cell_data,1))], 'FontSize',14);
    xlim([0 size(cell_data,2)+1]);ylim([min(median(X))-0.5 max(median(X))+2]);
    xticks(1:size(cell_data,2))
    xticklabels(markers);
    xtickangle(45)
    set(gca,'FontSize',8,'FontWeight','b');
    saveas(gcf,[nfolder '\barplots\PhenSig_' num2str(niche)],'png'); hold off;
    close
end
 
end
 
function [map, x_range, y_range] = kernel_density( x, y, bandwidth )
 
%Computes a gaussian kernel based density estimate of the points (x,y)
%Karthik Shekhar, MIT
 
 
%Params (can be changed)
width = 400; height = 400;
 
%Limits and deltax, deltay
limits(1) = min(x)-5;
limits(2) = max(x)+5;
limits(3) = min(y)-5;
limits(4) = max(y)+5;
 
deltax = (limits(2) - limits(1)) / width;
deltay = (limits(4) - limits(3)) / height;
 
if ~exist('bandwidth','var')
    disp('ERROR : must specify a value of the kernel bandwidth. Exiting...')
    return
end
 
map = zeros(height, width);
 
for ii = 0: height - 1
        yi = limits(3) + ii * deltay + deltay/2;
        for jj = 0 : width - 1
            xi = limits(1) + jj * deltax + deltax/2;
            
            
            dist2 = (x - xi).^2 + (y - yi).^2;
            dd = sum(exp(-dist2./(2*bandwidth^2))); 
 
            map(ii+1,jj+1) = (1/sqrt(2*pi*bandwidth^2)) * dd;
        end
end
            
    for ii=0:height-1
        y_range(ii+1) = limits(3) + ii * deltay + deltay/2;
    end
    for jj=0:width-1
        x_range(jj+1) = limits(1) + jj * deltax + deltax/2;
    end
    
    map = map/sum(map(:)); %Convert to PDF
    map(map < 5e-8)=0; 
end
 
function  centers=DetectPeaks(img, background,edg)
 
% Analyzes 2D images to detect peaks , finds x-y positions of peaks to 1 pixel accuracy
%
% Inputs:
%   img          The 2D data raw image 
%   background  A number between 0 and max(raw_image(:)) to remove  background
%   edg         A number>1 for skipping the first few and the last few 'edge' pixels
 
%
% Output:
%   centers        a 1xN vector of coordinates of peaks (x1,y1,x2,y2,...)
 
%This code is an adaptation of the skeleton http://www.mathworks.com/matlabcentral/fileexchange/37388-fast-2d-peak-finder/content/FastPeakFind.m
 
%% argcheck
 
if ~exist('img','var')
    disp('ERROR:Input image not provided. Exiting ... ');
    return
end
 
if isfloat(img) %For the case the input image is double, casting to uint16 keeps enough dynamic range while speeds up the code.
    img =  uint16( img.*2^16./(max(img(:))));
end
 
if (nargin < 2)
    background = (max([min(max(img,[],1))  min(max(img,[],2))])) ; 
end
 
if (nargin < 3)
    edg =3;
end
 
 
%% Analyze image
if any(img(:))   %for the case of non zero raw image
    
    
    % apply threshold
    img=img.*uint16(img>background);
    
    if any(img(:))    %for the case of the image is still non zero
        
      
        % img will be noisy on the edges, we'll skip 'edge' pixels.
        
        [x, y]=find(img(edg:size(img,1)-edg,edg:size(img,2)-edg));
        
        % initialize peak find outputs
        centers=[];
        x=x+edg-1;
        y=y+edg-1;
        for j=1:length(y)
            if (img(x(j),y(j))>=img(x(j)-1,y(j)-1 )) &&...
                    (img(x(j),y(j))>img(x(j)-1,y(j))) &&...
                    (img(x(j),y(j))>=img(x(j)-1,y(j)+1)) &&...
                    (img(x(j),y(j))>img(x(j),y(j)-1)) && ...
                    (img(x(j),y(j))>img(x(j),y(j)+1)) && ...
                    (img(x(j),y(j))>=img(x(j)+1,y(j)-1)) && ...
                    (img(x(j),y(j))>img(x(j)+1,y(j))) && ...
                    (img(x(j),y(j))>=img(x(j)+1,y(j)+1))
                
    
                centers = [centers ; x(j) ; y(j)];                    
            end
        end
        
        
    else % in case image after threshold is all zeros
        centers=[];
        return
    end
    
else % in case raw image is all zeros (dead event)
    centers=[];
    return
end
end
