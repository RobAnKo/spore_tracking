function [M,t] = directmatch(fquant1, fquant2,cfg,initial_t, beta)
%synopsis
%Mbin: binary correspondence matrix
%t: translational shift between frames

%cfg: configuration file
%initial_t: shift vector for first round of optimization
%beta: inverse temperature when doing softassign

var = cfg.EXP_VAR^2; %expected variance of cell positions of corresponding cells


%retrieve the cell positions
r1 = fquant1.spore.bf.center; %old positions
r1 = r1 - initial_t'; %correct with initial shift vector
r2 = fquant2.spore.bf.center; %new positions

%number of objects in resprective frames
f1max = size(r1,2);
f2max = size(r2,2);
           
%SOFTassign with high beta threshold: should be annealed when also
%determining A
if nargin<6
    beta = 1; 
end    
    
 
%initialize translation vector    
t = initial_t;
t_up = [0 0];
first_iter = true;
%update translation vector and correspondence matrix until shift is smaller than threshold
while norm(t_up) > 0.05 || first_iter
    Q = zeros(f1max,f2max);
    D = zeros(f1max,f2max);
    for i=1:f1max
        for j=1:f2max
              c=r1(:,i)-r2(:,j);
              D(i,j) = c'*c;
        end
    end
    
    % store initial distance matrix for later thresholding
    
    %if first_iter
    %    D_initial = D;
    %end
    
    %design of the Q-matrix for fixed A and t: scaled Euclidean distance with
    %bias towards matching if d1-d2 < var
    Q = -(D)/var+1;

    
    M = exp(beta*Q);

    %slack variable initialization for outlier row and column
    moutlier = 1/f1max * 0.01;
    m_outlier_row = ones(1,f2max) * moutlier; 
    m_outlier_col = ones(f1max,1) * moutlier; 
    [M, ~, ~] = normalize_m(M,m_outlier_col, m_outlier_row);
    
    %calculate the translation shift
    t_up(1) = 0;
    t_up(2) = 0;
    msum = 0;
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            t_up(1) = t_up(1)+M(i,j)*(r1(1,i)- r2(1,j));
            t_up(2) = t_up(2)+M(i,j)*(r1(2,i)- r2(2,j));
            msum = msum + M(i,j);
        end
    end

    if msum>0
        t_up=t_up./msum;
        t = t + t_up;
        
    end
    r1 = r1 - t_up';
    first_iter = false;
end

% last update of distance/correspondence matrix
Q = zeros(f1max,f2max);
D = zeros(f1max,f2max);
for i=1:f1max
    for j=1:f2max
          c=r1(:,i)-r2(:,j);
          D(i,j) = c'*c;
    end
end
Q = -(D)/var+1;
M = exp(beta*Q);
m_outlier_row = ones(1,f2max) * moutlier; 
m_outlier_col = ones(f1max,1) * moutlier; 
[M, m_outlier_col, m_outlier_row] = normalize_m(M,m_outlier_col, m_outlier_row);
for i = 1:f1max
    % check for high 'probability' of object having no corresponding object in other frame
    % OR objects where correspondence to more than one object is possible
    % and set to 0 -> no correspondence
    if m_outlier_col(i) > 0.05 || sum(M(i,:)) > 1.01 || max(M(i,:)) == 0
        M(i,:) = zeros(1,f2max);
    else
        % convert maximum values of each row/colums to 1, all others to 0
        M(i,:) = M(i,:) == max(M(i,:));
    end
end
for i = 1:f2max
    if m_outlier_row(i) > 0.05 || sum(M(:,i)) > 1.01 || max(M(:,i)) == 0
        M(:,i) = zeros(1,f1max);
    else
        M(:,i) = M(:,i) == max(M(:,i));
    end
end
              
              
