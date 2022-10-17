function [SHDT,CCstat,CCX,CCY] = shearletDescriptor(ST,Psi,L1,L2)
    % input: ST : shearlet coefficients computed by shearletTransformSpect
    %        Psi : spectrum of shearlets
    % output: SHDT : shearlet descriptor, 3D matrix
    %         CCstat : 
    
    % compute shearlet level and orientation
    [j,k] = arrayfun(@shearletScaleShear,1:size(ST,3));
    ST = ST(:,:,2:end);
    Psi = Psi(:,:,2:end);
    j = j(2:end);
    k = k(2:end);
    
    % find indices of shearlets in the 2 finest levels
    J = max(j);
    if J > 1
        DT_ind = find( j >= J-1 );
    else
        DT_ind = find( j == J );
    end
    use_ind = find( j >= J-2);
    DTj = j(DT_ind);
    
    % compute the weight of shearlets
    factor = sum(sum(abs(Psi).^2,1),2)/512^2;
    
    % normalized ST
    ST(:,:, use_ind ) = bsxfun(@times,ST(:,:,use_ind),...
                            reshape(2.^(j(use_ind)),[1,1,length(use_ind)])./sqrt(factor(use_ind)));
    
    % preload memory
    SHDT = zeros(size(ST,1),size(ST,2),length(DT_ind));
    
    %%
    % go through each shearlet and find the local max response
    % at each pixel in the shearlet parameter space
    DTi = 1;
    ST = abs(ST);

% %     gaussf = fspecial('gaussian',10,1);
% %     gaussf = imresize(gaussf,[10,20]);
% %     for i = 1:8
% %         ST(:,:,DT_ind(i)) = imfilter(ST(:,:,DT_ind(i)),imrotate(gaussf,22.5*(i-1),'bicubic'));
% % %         imshow(ST(:,:,DT_ind(i)),[]);pause;
% %     end
% %     gaussf = imresize(gaussf,[10,10]);
% %     for i = 9:24
% %         ST(:,:,DT_ind(i)) = imfilter(ST(:,:,DT_ind(i)),gaussf);
% % %         imshow(ST(:,:,DT_ind(i)),[]);pause;
% %     end
% %     
% %     medST = arrayfun(@(x)medfilt2(ST(:,:,x),[1,1]*2^(J-j(x))+1),1:size(ST,3),'UniformOutput',0);
% %     medST = cellfun(@(x)(x+mean(x(:))),medST,'UniformOutput',0);
% %     medST = cat(3,medST{:});
% %     ST = ST - medST*.9;
% %     ST = ST.*(ST>0);
% %     ST = ST./repmat(sum(ST,3)+1000,[1,1,size(ST,3)]);
% % 
% %     figure;for i = 1:24, imshow(ST(:,:,4+i),[]);pause;end
    
    % focus on the central part of the image
    if ~exist('L1','var')
        L1 = floor(size(ST,1))/4;
    end
    if ~exist('L2','var')
        L2 = L1;
    end
    L = floor(size(ST,1)/2);
    min_resp = max(reshape(ST(L-L1:L+L1,L-L2:L+L2,use_ind),[],1))/3;
    for ind_j = 0:1
        curr_j = find(j == J + ind_j -1);
        parent_j = find(j == J + ind_j - 2);
        parent_j = [parent_j;parent_j(2:end) parent_j(1)];
        parent_j = parent_j(:)';
        if ~ind_j
            child_j = find( j == J + ind_j);
            child_j = reshape([child_j(end) child_j(1:end-1)],2,[]);
        else
            child_j = [];
        end
        neigh_j = [curr_j(end) curr_j(1:end-1); curr_j(2:end) curr_j(1)];
        for i = 1:length(curr_j)
            if ~ind_j
                subST = ST(:,:,[curr_j(i),neigh_j(:,i)',parent_j(i),child_j(:,i)']);
            else
                subST = ST(:,:,[curr_j(i),neigh_j(:,i)',parent_j(i)]);
            end
            [~,max_ind] = max(subST,[],3);
            SHDT(:,:,DTi) = (max_ind <= 3) & (ST(:,:,curr_j(i)) > min_resp);
            DTi = DTi + 1;
        end
    end
    
    %%
    % identify large connected components, w.r.t. the size of the image
    con_cmp_size = size(ST,1)/20; % at least one-third of the size of the image
%     % weak-connect ==> strong-connect
%     neigh_conn = conndef(2,'minimal');
%     SHDT = arrayfun(@(x)imfilter(SHDT(:,:,x),neigh_conn),1:size(SHDT,3),...
%                     'UniformOutput',0);
    SHDT = mat2cell(SHDT,size(SHDT,1),size(SHDT,2),ones(size(SHDT,3),1));
%     for i = 1:24;imshow(label2rgb(bwlabel(SHDT{i},8)),[]);pause;end
    % find connected components
    CC = cellfun(@(x)bwconncomp(x,4),SHDT,'UniformOutput',0);
    % exclude small components, sort components in decreasing size
    for i = 1:length(CC)
        CCsize = cellfun(@length,CC{i}.PixelIdxList);
        CC{i}.PixelIdxList = CC{i}.PixelIdxList(CCsize > con_cmp_size);
        CC{i}.NumObjects = length(CC{i}.PixelIdxList);
%         figure;imshow(SHDT{i});
%         pause;
    end
    % extract statistics of each component
    CCstat = cellfun(@(x)regionprops(x,'Centroid','Orientation','MajorAxisLength'),...
                        CC,'UniformOutput',0);
    
    % compute line segments of each component
    CCX = [];
    CCY = [];
    CCorient = [];% orientation
    CCr = [];% radius
    CCj = [];% scale
    for i = 1:length(CC)
        CClength = [CCstat{i}(:).MajorAxisLength];
        % discard components not elongate enough
        [CClength,CCind] = sort(CClength,'descend');
        CC{i}.NumObjects = find(CClength >= con_cmp_size,1,'last');
        if isempty(CC{i}.NumObjects)
            CC{i}.NumObjects = 0;
            continue;
        end
        CCind = CCind(1:CC{i}.NumObjects);
        CCstat{i} = CCstat{i}(CCind);
        CC{i}.PixelIdxList = CC{i}.PixelIdxList(CCind);
        CClength = CClength(1:CC{i}.NumObjects);
        % compute angle and the center of line segment
        CCangle = -[CCstat{i}(:).Orientation]/180*pi;
%         normal_vector = median([-sin(CCangle); cos(CCangle)],2);
%         tang_vector = [normal_vector(2); -normal_vector(1)];
        CCloc = cat(1,CCstat{i}.Centroid);
        if isempty(CCloc)
            continue;
        end
%         % discard redundant components
%         pt = 1;
%         while pt <= length(CC{i}.PixelIdxList)
%             % compute distance between centroids of components
%             dist = bsxfun(@minus,CCloc(pt+1:end,:),CCloc(pt,:));
%             n_dist = abs(dist*normal_vector);
%             t_dist = abs(dist*tang_vector);
%             % find redundant components
%             repind = find((n_dist < 2^(J-DTj(i)+1)*2 & t_dist < 2^(2*(J-DTj(i)+1)))|...
%                 (n_dist < 2^(J-DTj(i)+1) & t_dist < 2^(2*(J-DTj(i)+1))*2)) + pt;
%             CC{i}.PixelIdxList(repind) = [];
%             CCloc(repind,:) = [];
%             CClength(repind) = [];
%             pt = pt + 1;
%         end
        CC{i}.NumObjects = length(CC{i}.PixelIdxList);
        CCangle = regionprops(CC{i},'Orientation');
        CCangle = -[CCangle(:).Orientation]/180*pi;
        dX = CClength.*cos(CCangle)/2;
        X = [ CCloc(:,1)' + dX; CCloc(:,1)' - dX];
        dY = CClength.*sin(CCangle)/2;
        Y = [ CCloc(:,2)' + dY; CCloc(:,2)' - dY];
% %         figure;imshow(SHDT{i});hold on;
% %         line(X,Y,'LineWidth',2,'Color','r');hold off
% %         pause;
        CCX = [CCX X];
        CCY = [CCY Y];
        CCorient = [CCorient CCangle];
        CCr = [CCr CClength];
        CCj = [CCj DTj(i)*ones(1,size(X,2))];
    end

    % discard redundant componets: close in orientation and distance
    % compute the statistics of each component
    [CCr,CCind] = sort(CCr,'descend');
    CCX = CCX(:,CCind);
    CCY = CCY(:,CCind);
    CCorient = CCorient(:,CCind);
    CCr = CCr(:,CCind);
    CCloc = arrayfun(@(x)cat(1,CCstat{x}.Centroid), 1:length(CCstat),'UniformOutput',0);
    CCloc = cat(1,CCloc{:});
    CCloc = CCloc(CCind,:);
    CCj = CCj(CCind);
    
    pt = 1;
    length(CCr)
    while pt < length(CCr)
            dist = bsxfun(@minus,CCloc(pt+1:end,:),CCloc(pt,:));
            normal_vector = [-sin(CCorient(pt)); cos(CCorient(pt))];
            tang_vector = [normal_vector(2); -normal_vector(1)];
           n_dist = abs(dist*normal_vector);
            t_dist = abs(dist*tang_vector);
            % find redundant components
            repind = find(n_dist < 2^(J-CCj(pt)+1)*2 & t_dist < CCr(pt) & cos(CCorient(pt+1:end)'-CCorient(pt)) > .9) + pt;
            CCloc(repind,:) = [];
            CCr(repind) = [];
            CCorient(repind) = [];
            CCX(:,repind) = [];
            CCY(:,repind) = [];
            pt = pt + 1;
    end
    length(CCr)
    
    flag = cellfun(@(x)(x.NumObjects > 0),CC); 
    
    CCstat = cat(1,CCstat{flag});
    
    figure;%imshow(ones(size(SHDT{1})),[]);hold on;
    imshow(sum(cat(3,SHDT{:}),3),[0,2]);hold on
    line(CCX,CCY,'LineWidth',2,'Color','r');
    hold off