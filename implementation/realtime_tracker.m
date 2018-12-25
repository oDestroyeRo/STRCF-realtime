function params = realtime_tracker(params, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get params.sequence info
%[params.seq, im] = get_params.sequence_info(params.seq);
%params.seq = params.seq;
if ~(params.seq.isFirst)
    params.seq.init_sz = [params.seq.init_rect(1,4), params.seq.init_rect(1,3)];
    params.seq.init_pos = [params.seq.init_rect(1,2), params.seq.init_rect(1,1)] + (params.seq.init_sz - 1)/2;
    params.im = params.seq.img;
    %params = rmfield(params, 'params.seq');
    if isempty(params.im)
        %params.seq.rect_position = [];
        %[~, results] = get_params.sequence_results(params.seq);
        return;
    end

    % Init position
    params.pos = params.seq.init_pos(:)';
    params.target_sz = params.seq.init_sz(:)';
    params.init_sz = params.target_sz;

    % Feature settings
    params.features = params.t_features;

    % Set default parameters
    params = init_default_params(params);
    % Global feature parameters
    if isfield(params, 't_global')
        params.global_fparams = params.t_global;
    else
        params.global_fparams = [];
    end

    params.global_fparams.use_gpu = params.use_gpu;
    params.global_fparams.gpu_id = params.gpu_id;

    % Define data types
    if params.use_gpu
        params.data_type = zeros(1, 'single', 'gpuArray');
    else
        params.data_type = zeros(1, 'single');
    end
    params.data_type_complex = complex(params.data_type);

    params.global_fparams.data_type = params.data_type;

    % Load learning parameters
    params.admm_max_iterations = params.max_iterations;
    params.init_penalty_factor = params.init_penalty_factor;
    params.max_penalty_factor = params.max_penalty_factor;
    params.penalty_scale_step = params.penalty_scale_step;
    params.temporal_regularization_factor = params.temporal_regularization_factor; 

    params.init_target_sz = params.target_sz;

    % Check if color image
    if size(params.im,3) == 3
        if all(all(params.im(:,:,1) == params.im(:,:,2)))
            params.is_color_image = false;
        else
            params.is_color_image = true;
        end
    else
        params.is_color_image = false;
    end

    if size(params.im,3) > 1 && params.is_color_image == false
        params.im = params.im(:,:,1);
    end

    % Check if mexResize is available and show warning otherwise.
    params.use_mexResize = true;
    params.global_fparams.use_mexResize = true;
    try
        [~] = mexResize(ones(5,5,3,'uint8'), [3 3], 'auto');
    catch
        params.use_mexResize = false;
        params.global_fparams.use_mexResize = false;
    end

    % Calculate search area and initial scale factor
    params.search_area = prod(params.init_target_sz * params.search_area_scale);
    if params.search_area > params.max_image_sample_size
        params.currentScaleFactor = sqrt(params.search_area / params.max_image_sample_size);
    elseif params.search_area < params.min_image_sample_size
        params.currentScaleFactor = sqrt(params.search_area / params.min_image_sample_size);
    else
        params.currentScaleFactor = 1.0;
    end

    % target size at the initial scale
    params.base_target_sz = params.target_sz / params.currentScaleFactor;

    % window size, taking padding into account
    switch params.search_area_shape
        case 'proportional'
            params.img_sample_sz = floor(params.base_target_sz * params.search_area_scale);     % proportional area, same aspect ratio as the target
        case 'square'
            params.img_sample_sz = repmat(sqrt(prod(params.base_target_sz * params.search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
        case 'fix_padding'
            params.img_sample_sz = params.base_target_sz + sqrt(prod(params.base_target_sz * params.search_area_scale) + (params.base_target_sz(1) - params.base_target_sz(2))/4) - sum(params.base_target_sz)/2; % const padding
        case 'custom'
            params.img_sample_sz = [params.base_target_sz(1)*2 params.base_target_sz(2)*2];
    end

    [params.features, params.global_fparams, params.feature_info] = init_features(params.features, params.global_fparams, params.is_color_image, params.img_sample_sz, 'exact');

    % Set feature info
    params.img_support_sz = params.feature_info.img_support_sz;
    params.feature_sz = unique(params.feature_info.data_sz, 'rows', 'stable');
    params.feature_cell_sz = unique(params.feature_info.min_cell_size, 'rows', 'stable');
    params.num_feature_blocks = size(params.feature_sz, 1);

    % Get feature specific parameters
    params.feature_extract_info = get_feature_extract_info(params.features);

    % Size of the extracted feature maps
    params.feature_sz_cell = mat2cell(params.feature_sz, ones(1,params.num_feature_blocks), 2);
    params.filter_sz = params.feature_sz;
    params.filter_sz_cell = permute(mat2cell(params.filter_sz, ones(1,params.num_feature_blocks), 2), [2 3 1]);

    % The size of the label function DFT. Equal to the maximum filter size
    [params.output_sz, params.k1] = max(params.filter_sz, [], 1);
    params.k1 = params.k1(1);

    % Get the remaining block indices
    params.block_inds = 1:params.num_feature_blocks;
    params.block_inds(params.k1) = [];

    % Construct the Gaussian label function
    params.yf = cell(numel(params.num_feature_blocks), 1);
    for i = 1:params.num_feature_blocks
        params.sz = params.filter_sz_cell{i};
        params.output_sigma = sqrt(prod(floor(params.base_target_sz/params.feature_cell_sz(i)))) * params.output_sigma_factor;
        params.rg           = circshift(-floor((params.sz(1)-1)/2):ceil((params.sz(1)-1)/2), [0 -floor((params.sz(1)-1)/2)]);
        params.cg           = circshift(-floor((params.sz(2)-1)/2):ceil((params.sz(2)-1)/2), [0 -floor((params.sz(2)-1)/2)]);
        [params.rs, params.cs]     = ndgrid(params.rg,params.cg);
        params.y            = exp(-0.5 * (((params.rs.^2 + params.cs.^2) / params.output_sigma^2)));
        params.yf{i}           = fft2(params.y); 
    end

    % Compute the cosine windows
    params.cos_window = cellfun(@(sz) hann(params.sz(1))*hann(params.sz(2))', params.feature_sz_cell, 'uniformoutput', false);

    % Define spatial regularization windows
    params.reg_window = cell(params.num_feature_blocks, 1);
    for i = 1:params.num_feature_blocks
        params.reg_scale = floor(params.base_target_sz/params.feature_downsample_ratio(i));
        params.use_sz = params.filter_sz_cell{i};    
        params.reg_window{i} = ones(params.use_sz) * params.reg_window_max;
        params.range = zeros(numel(params.reg_scale), 2);

        % determine the target center and range in the regularization windows
        for j = 1:numel(params.reg_scale)
            params.range(j,:) = [0, params.reg_scale(j) - 1] - floor(params.reg_scale(j) / 2);
        end
        params.center = floor((params.use_sz + 1)/ 2) + mod(params.use_sz + 1,2);
        params.range_h = (params.center(1)+ params.range(1,1)) : (params.center(1) + params.range(1,2));
        params.range_w = (params.center(2)+ params.range(2,1)) : (params.center(2) + params.range(2,2));

        params.reg_window{i}(params.range_h, params.range_w) = params.reg_window_min;
    end

    % Pre-computes the grid that is used for socre optimization
    params.ky = circshift(-floor((params.filter_sz_cell{1}(1) - 1)/2) : ceil((params.filter_sz_cell{1}(1) - 1)/2), [1, -floor((params.filter_sz_cell{1}(1) - 1)/2)]);
    params.kx = circshift(-floor((params.filter_sz_cell{1}(2) - 1)/2) : ceil((params.filter_sz_cell{1}(2) - 1)/2), [1, -floor((params.filter_sz_cell{1}(2) - 1)/2)])';
    params.newton_iterations = params.newton_iterations;

    % Use the translation filter to estimate the scale
    params.nScales = params.number_of_scales;
    params.scale_step = params.scale_step;
    params.scale_exp = (-floor((params.nScales-1)/2):ceil((params.nScales-1)/2));
    params.scaleFactors = params.scale_step .^ params.scale_exp;

    if params.nScales > 0
        %force reasonable scale changes
        params.min_scale_factor = params.scale_step ^ ceil(log(max(5 ./ params.img_support_sz)) / log(params.scale_step));
        params.max_scale_factor = params.scale_step ^ floor(log(min([size(params.im,1) size(params.im,2)] ./ params.base_target_sz)) / log(params.scale_step));
    end
    
    % Define the learning variables
    params.f_pre_f = cell(params.num_feature_blocks, 1);
    params.cf_f = cell(params.num_feature_blocks, 1);

    % Allocate
    params.scores_fs_feat = cell(1,1,params.num_feature_blocks);
    
    params.seq.isFirst = true;
end

if (params.seq.isFirst)
    % Read image
    if params.seq.isOldPos == true
        %[params.seq, im] = get_params.sequence_frame(params.seq);
        %params.im = snapshot(handles.cam);
        %params.im = flip(params.im,2);
        %params.seq.img_real = params.im;
        %params.im = rgb2gray(params.im);
        if isempty(params.im)
           return
        end
        if size(params.im,3) > 1 && params.is_color_image == false
            params.im = params.im(:,:,1);
        end
    else
        params.seq.isInit = true;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Target localization step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Do not estimate translation and scaling on the first frame, since we 
    % just want to initialize the tracker there
    if params.seq.isOldPos == true
        params.old_pos = inf(size(params.pos));
        params.iter = 1;
        
        %translation search
        while params.iter <= params.refinement_iterations && any(params.old_pos ~= params.pos)
            % Extract features at multiple resolutions
            params.sample_pos = round(params.pos);
            params.sample_scale = params.currentScaleFactor*params.scaleFactors;
            params.xt = extract_features(params.im, params.sample_pos, params.sample_scale, params.features, params.global_fparams, params.feature_extract_info);
                                    
            % Do windowing of features
            
            params.xtw = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), params.xt, params.cos_window, 'uniformoutput', false);
            
            % Compute the fourier series
            params.xtf = cellfun(@fft2, params.xtw, 'uniformoutput', false);
                        
            % Compute convolution for each feature block in the Fourier domain
            % and the sum over all blocks.
            params.scores_fs_feat{params.k1} = gather(sum(bsxfun(@times, conj(params.cf_f{params.k1}), params.xtf{params.k1}), 3));
            params.scores_fs_sum = params.scores_fs_feat{params.k1};
            for k = params.block_inds
                params.scores_fs_feat{k} = gather(sum(bsxfun(@times, params.conj(params.cf_f{k}), params.xtf{k}), 3));
                params.scores_fs_feat{k} = resizeDFT2(params.scores_fs_feat{k}, params.output_sz);
                params.scores_fs_sum = params.scores_fs_sum +  params.scores_fs_feat{k};               
            end
            % Also sum over all feature blocks.
            % Gives the fourier coefficients of the convolution response.
            params.scores_fs = permute(gather(params.scores_fs_sum), [1 2 4 3]);
            
            params.responsef_padded = resizeDFT2(params.scores_fs, params.output_sz);
            params.response = ifft2(params.responsef_padded, 'symmetric');
            [params.disp_row, params.disp_col, params.sind] = resp_newton(params.response, params.responsef_padded, params.newton_iterations, params.ky, params.kx, params.output_sz);
                        
            % Compute the translation vector in pixel-coordinates and round
            % to the closest integer pixel.
            params.translation_vec = [params.disp_row, params.disp_col] .* (params.img_support_sz./params.output_sz) * params.currentScaleFactor * params.scaleFactors(params.sind);            
            params.scale_change_factor = params.scaleFactors(params.sind);
            
            % update position
            params.old_pos = params.pos;
            params.pos = params.sample_pos + params.translation_vec;
            if params.clamp_position
                params.pos = max([1 1], min([size(params.im,1) size(params.im,2)], params.pos));
            end
                        
            % Update the scale
            params.currentScaleFactor = params.currentScaleFactor * params.scale_change_factor;
            
            % Adjust to make sure we are not to large or to small
            if params.currentScaleFactor < params.min_scale_factor
                params.currentScaleFactor = params.min_scale_factor;
            elseif params.currentScaleFactor > params.max_scale_factor
                params.currentScaleFactor = params.max_scale_factor;
            end
            
            params.iter = params.iter + 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Model update step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract image region for training sample
    params.sample_pos = round(params.pos);

    params.xl = extract_features(params.im, params.sample_pos, params.currentScaleFactor, params.features, params.global_fparams, params.feature_extract_info);
    
    % do windowing of features
    params.xlw = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), params.xl, params.cos_window, 'uniformoutput', false);
    % compute the fourier series
    params.xlf = cellfun(@fft2, params.xlw, 'uniformoutput', false);
    % train the CF model for each feature
    for k = 1: numel(params.xlf)
        params.model_xf = params.xlf{k};

        if (params.seq.isInit == true)
            params.f_pre_f{k} = zeros(size(params.model_xf));
            params.mu = 0;
            params.seq.isInit = false;
        else
            params.mu = params.temporal_regularization_factor(k);
        end
        
        % intialize the variables
        params.f_f = single(zeros(size(params.model_xf)));
        params.g_f = params.f_f;
        params.h_f = params.f_f;
        params.gamma  = params.init_penalty_factor(k);
        params.gamma_max = params.max_penalty_factor(k);
        params.gamma_scale_step = params.penalty_scale_step(k);
        
        % use the GPU mode
        if params.use_gpu
            params.model_xf = gpuArray(params.model_xf);
            params.f_f = gpuArray(params.f_f);
            params.f_pre_f{k} = gpuArray(params.f_pre_f{k});
            params.g_f = gpuArray(params.g_f);
            params.h_f = gpuArray(params.h_f);
            params.reg_window{k} = gpuArray(params.reg_window{k});
            params.yf{k} = gpuArray(params.yf{k});
        end

        % pre-compute the variables
        params.T = prod(params.output_sz);
        params.S_xx = sum(conj(params.model_xf) .* params.model_xf, 3);
        params.Sf_pre_f = sum(conj(params.model_xf) .* params.f_pre_f{k}, 3);
        params.Sfx_pre_f = bsxfun(@times, params.model_xf, params.Sf_pre_f);

        % solve via ADMM algorithm
        params.iter = 1;
        while (params.iter <= params.admm_max_iterations)

            % subproblem f
            params.B = params.S_xx + params.T * (params.gamma + params.mu);
            params.Sgx_f = sum(conj(params.model_xf) .* params.g_f, 3);
            params.Shx_f = sum(conj(params.model_xf) .* params.h_f, 3);
 
            params.f_f = ((1/(params.T*(params.gamma + params.mu)) * bsxfun(@times,  params.yf{k}, params.model_xf)) - ((1/(params.gamma + params.mu)) * params.h_f) +(params.gamma/(params.gamma + params.mu)) * params.g_f) + (params.mu/(params.gamma + params.mu)) * params.f_pre_f{k} - ...
                bsxfun(@rdivide,(1/(params.T*(params.gamma + params.mu)) * bsxfun(@times, params.model_xf, (params.S_xx .*  params.yf{k})) + (params.mu/(params.gamma + params.mu)) * params.Sfx_pre_f - ...
                (1/(params.gamma + params.mu))* (bsxfun(@times, params.model_xf, params.Shx_f)) +(params.gamma/(params.gamma + params.mu))* (bsxfun(@times, params.model_xf, params.Sgx_f))), params.B);

            %   subproblem g
            params.g_f = fft2(argmin_g(params.reg_window{k}, params.gamma, real(ifft2(params.gamma * params.f_f+ params.h_f)), params.g_f));

            %   update h
            params.h_f = params.h_f + (params.gamma * (params.f_f - params.g_f));

            %   update gamma
            params.gamma = min(params.gamma_scale_step * params.gamma, params.gamma_max);
            
            params.iter = params.iter+1;
        end
        
        % save the trained filters
        params.f_pre_f{k} = params.f_f;
        params.cf_f{k} = params.f_f;
    end  
            
    % Update the target size (only used for computing output box)
    params.target_sz = params.base_target_sz * params.currentScaleFactor;
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Visualization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % visualization
    if params.visualization
        params.rect_position_vis = [params.pos([2,1]) - (params.target_sz([2,1]) - 1)/2, params.target_sz([2,1])];
    end
    params.seq.isFirst = true;
    params.seq.isOldPos = true;
end

