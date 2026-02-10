function [HMRFseg, settings] = HMRFseg4img(img, nucleus_mask, nclust, beta0, mineps)
    % ----------------------- parameters settings -------------------------
    img = double(img);
    nucleus_mask = logical(nucleus_mask); mask = nucleus_mask(:);
    maxiter = 200;
%     mineps = 10^(-7);
    varfixed = true;
    inforce = true;
    
    % ---------- Hidden Markov Random Field (HMRF) model --------------
    % initiation
    X = size(nucleus_mask, 1); Y = size(nucleus_mask, 2);
    normf = quantile(img(mask), 0.9995);
    img = img(:)/normf; 
    img_class = zeros(size(mask)); % from 1 to nclust
    qu = quantile(img(mask), (1:(nclust-1))/nclust);
    for i = 1:(nclust-1)
        img_class((img>qu(i))&mask) = i;
    end
    img_class(mask) = img_class(mask) + 1;
    beta = diag(ones(1, nclust)*beta0);
    
    mu = zeros(nclust, 1);
    sigma = ones(nclust, 1);
    
    for i = 1:nclust
        mu(i) = mean(img(img_class==i));
        sigma(i) = std(img(img_class==i));
    end
    if varfixed
        sigma = ones(nclust, 1);
        sigma = std(img(mask)-mu(img_class(mask)))*sigma;
    end
    loglik = zeros(nclust, 1);
    counter = 0;
    criterium = true;
    
    % interation
    while criterium
        counter = counter + 1;
        disp(['Iteration: ', num2str(counter)]);
    
        for xidx = 1:X
            for yidx = 1:Y
                idx = (yidx-1)*X + xidx;
                if mask(idx) 
                    for i = 1:nclust
                        loglik(i) = -0.5*(img(idx)-mu(i))^2/sigma(i)/sigma(i);
                    end
                    if ~(xidx == 1)
                        nid = (yidx-1)*X+xidx-1;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    if ~(xidx == X)
                        nid = (yidx-1)*X+xidx+1;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    if ~(yidx == 1)
                        nid = (yidx-2)*X+xidx;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    if ~(yidx == Y)
                        nid = yidx*X+xidx;
                        if mask(nid)
                        for i = 1:nclust
                            loglik(i) = loglik(i) + mask(nid)*beta(img_class(nid), i);
                        end
                        end
                    end
                    [~, temp] = max(loglik);
                    img_class(idx) = temp;
                end
            end
        end
    
        temp_nclust = nclust;
        for i = nclust:-1:1
            if sum(img_class == i) == 0
                disp(["class ", num2str(i), " removed."]);
                img_class(img_class>i) = img_class(img_class>i)-1;
                temp_nclust = temp_nclust - 1;
            end
        end
        oldmu = mu(1:temp_nclust);
        mu = zeros(temp_nclust, 1);
        sigma = ones(temp_nclust, 1);
        for i = 1:temp_nclust
            mu(i) = mean(img(img_class==i));
            sigma(i) = std(img(img_class==i));
        end
        if varfixed
            sigma = ones(nclust, 1);
            sigma = std(img(mask)-mu(img_class(mask)))*sigma;
        end
        % re-arange img_class based on intensity
        if ~issorted(mu)
            [~, mu_idx] = sort(mu); 
            mu = mu(mu_idx); sigma = sigma(mu_idx);
            for i = 1:length(mu_idx)
                if mu(i) ~= i
                    img_class(img_class==mu(i)) = i;
                end
            end
        end
        [~,~,ic] = unique(img_class);
        a_counts = accumarray(ic,1);
        sigma(sigma==0) = 1e-6;
        disp("Class number: ");
        disp(num2str(a_counts));
        disp("Class mean:");
        disp(num2str(mu));
        disp("Class sigma:");
        disp(num2str(sigma));
    
        if (counter == maxiter)||(sum((mu-oldmu).^2)<mineps)
            criterium = false;
        end
    
        if inforce
            while temp_nclust < nclust
                criterium = true;
                disp(["inforce nclust: ", num2str(nclust)]);
                if ~varfixed
                    [~, w] = max(sigma);
                else
                    [~,~,ic] = unique(img_class);
                    a_counts = accumarray(ic,1);
                    [~, w] = max(a_counts(2:end));
                end
                img_class(img_class>w) = img_class(img_class>w) + 1;
                nn = sum(img_class == w);
                disp(['Class with max sigma has pixel number: ', num2str(nn)]);
                img_class((img_class==w)&(img>mu(w))) = img_class((img_class==w)&(img>mu(w))) + 1;
                temp_nclust = temp_nclust + 1;
                mu = zeros(temp_nclust, 1);
                sigma = ones(temp_nclust, 1);
                for i = 1:temp_nclust
                    mu(i) = mean(img(img_class==i));
                    sigma(i) = std(img(img_class==i));
                end
                if varfixed
                    sigma = ones(nclust, 1);
                    sigma = std(img(mask)-mu(img_class(mask)))*sigma;
                end
                sigma(sigma==0) = 1e-6;
            end
        end
    end
    disp(['Stop at iteration: ', num2str(counter)]);

    [~,~,ic] = unique(img_class);
    a_counts = accumarray(ic,1); a_counts = a_counts((end-nclust+1):end);

    HMRFseg.img = reshape(img, size(nucleus_mask))*normf;
    HMRFseg.img_class = reshape(img_class, size(nucleus_mask));
    HMRFseg.nucleus_mask = nucleus_mask; % HMRFseg.mask = mask; 
    HMRFseg.mu = mu*normf;
    HMRFseg.sigma = sigma*normf;
    HMRFseg.a_counts = a_counts;

    settings.maxiter = maxiter; settings.normf = normf;
    settings.mineps = mineps; settings.nclust = nclust; settings.beta0 = beta0;
    settings.varfixed = varfixed; settings.inforce = inforce; 
end