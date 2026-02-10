function [img_deconv, img_reconv] = deconv_fixed_iter(img, psf, pad, iter)

    img_pad = padarray(img, [pad pad], 'symmetric');

    disp(['Iteration: ', num2str(iter)])
    
    % RL deconvolution
    Jpad = deconvlucy(img_pad, psf, iter);

    % Reblur
    Reblurred = imfilter(Jpad, psf, "symmetric", "conv");

    idx = repmat({':'}, 1, ndims(Jpad));
    idx{1} = pad+1 : size(Jpad,1)-pad;
    idx{2} = pad+1 : size(Jpad,2)-pad;

    % Crop
    img_deconv = Jpad(idx{:});
    img_reconv = Reblurred(idx{:});

end