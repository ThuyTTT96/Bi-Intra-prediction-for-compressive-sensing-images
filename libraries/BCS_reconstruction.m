function [reconstructed_image, MIN_ENERGY, L1_NORM] = BCS_reconstruction(y_deresidual, theta, phi, image_transformation_algorithm, image_reconstruction_algorithm, sub_pixels, opts)
    %___Initial guess = min energy IDCT___
    MIN_ENERGY = theta'*(y_deresidual);

    %___MINIMIZATION___
    %___SOLVE THE LP (L1-NORM SOLUTION)___
    switch image_reconstruction_algorithm
        case 'l1_eq_pd_spgl'
            L1_NORM = spg_bp(theta, ...
                             y_deresidual, ...
                             opts);
        case 'l1_eq_pd'
            L1_NORM = l1eq_pd(MIN_ENERGY, ...
                              theta, ...
                              [], ...
                              y_deresidual, ...
                              1e-3); % L1-magic toolbox
        case 'SPGL1'
            opts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
            L1_NORM = spg_bp(phi, y_deresidual, opts);
        case 'l1_qc_logbarrier'
            % number of observations to make
            sigma = 0.005;
            K = 1;
            epsilon = sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K));
            L1_NORM = l1qc_logbarrier(MIN_ENERGY, theta, [], y_deresidual, epsilon, 1e-3);
        case 'D-AMP' 
            %with denoise 
            
    end
    
    %___L1-NORM IMAGE RECONSTRUCTIONS with IFWHT WAVELET___
    switch image_transformation_algorithm
        case 'idct'
            final_image_transformation = idct2(L1_NORM);
        case 'ifwht'
            final_image_transformation = ifwht(L1_NORM);
        case 'idwt'
            % X contains the loaded image. 
            sX = size(L1_NORM);
            % Perform single-level decomposition 
            % of X using db4. 
            [cA1,cH1,cV1,cD1] = dwt2(L1_NORM,'db4');
            % Invert directly decomposition of X 
            % using coefficients at level 1. 
            final_image_transformation = idwt2(cA1,cH1,cV1,cD1,'db4',sX);
        otherwise
            disp('contact Jay')
    end
    
    %MIN_ENERGY = reshape(MIN_ENERGY.',sub_pixels,[])';
    %L1_NORM = reshape(L1_NORM.',sub_pixels,[])';
    reconstructed_image = reshape(final_image_transformation.',sub_pixels,[])'; %raster
    %reconstructed_image = izigzag(final_image_transformation, 4, 4)';
end