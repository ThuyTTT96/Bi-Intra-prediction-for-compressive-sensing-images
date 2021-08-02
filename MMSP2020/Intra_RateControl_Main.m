clear;
close all;
clc;
profile on
UpFolder = fileparts(pwd);
addpath(fullfile(UpFolder, 'libraries/l1magic/Optimization'));
addpath(fullfile(UpFolder, 'libraries/l1magic/Measurements'));
addpath(fullfile(UpFolder, 'libraries/l1magic/Data'));
addpath(fullfile(UpFolder, 'libraries/spgl1-1.9'));
addpath(fullfile(UpFolder, 'libraries'));
addpath(fullfile(UpFolder,'results'));
addpath(genpath('C:/Users/Thuy Tran/Documents/Research/Frame/WILD')); %replace some personal paths like this
addpath(genpath('C:/Users/Thuy Tran/Documents/Research/FPGA_ASIC/gray_compressed_sensing_walsh_ifwht_intra_prediction/y_reconstruct'));
 
%camList = webcamlist;
%Connect to the webcam.
%cam = webcam(1);
%preview(cam);

opts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output

imageOriginalPath = 'C:/Users/Thuy Tran/Documents/Research/Frame/WILD';
imageFiles = [dir(fullfile(imageOriginalPath,'*png'));
              dir(fullfile(imageOriginalPath,'*tiff'));
              dir(fullfile(imageOriginalPath,'*jpg'));
              dir(fullfile(imageOriginalPath,'*bmp'));
              dir(fullfile(imageOriginalPath,'*pgm'))
              dir(fullfile(imageOriginalPath,'*mat'))];
numFiles = length(imageFiles)
%___SIMULATION SETUPS___
sub_pixels       = 16;
n                = sub_pixels*sub_pixels; % NOTE: small error still present after increasing m;
bpp_buffer       = 0;

measurement_matrix_lists        = [128];
measurement_matrix_construction = 'binary_walsh_hadamard';
image_reconstruction_algorithm  = 'l1_eq_pd';
image_transformation_algorithm  = 'ifwht';
color_mode                      = 'gray';
quartus_interface               = 'off ';

stack_size = 3;

for matrix_loop = 1:length(measurement_matrix_lists)
    reset          = 1;
    switch measurement_matrix_lists(matrix_loop)
        case {256}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 1;
        case {192}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.75;
        case {128}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.50;
        case {64}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.25;
        case {32}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.125;
        case {16}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.0625;
        case {8}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.5;
        case {4}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.015625;
        case {2}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.0078125;
    end

        switch measurement_matrix_construction
            case 'binary_random'
                phi                      = randi([0, 1], [m, n]); % This will give you a friendly measurement matrix (M must smaller than N)
            case 'binary_walsh_hadamard'
                hadamard_matrix          = hadamard(n);
                HadIdx                   = 0:n-1;                         % Hadamard index
                M                        = log2(n)+1;                     % Number of bits to represent the index
                binHadIdx                = fliplr(dec2bin(HadIdx,M))-'0'; % Bit reversing of the binary index
                binSeqIdx                = zeros(n,M-1);                  % Pre-allocate memory
                for k = M:-1:2
                    % Binary sequency index
                    binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
                end
                SeqIdx                   = binSeqIdx*pow2((M-1:-1:0)');    % Binary to integer sequency index
                walshMatrix              = hadamard_matrix(SeqIdx+1,:); % 1-based indexing
                phi                      = max(walshMatrix(1:m,1:n), 0);
                      
            case 'guassian'
                phi                      = randn(m,n);
            case 'binary_checkerboard'
                n_n = sub_pixels*sub_pixels;    % the size of the matrix
                i   = 1 : n_n;                  % the index list
                p   = mod(i, 2);                % calculate parity
                r   = repmat(p,[n_n 1]);        % expand the parities to rows
                c   = repmat(p',[1 n_n]);       % expand the parities to columns
                yay = xor(r, c);                % apply the xor operator
                phi = double(yay(1:m,1:n));
            case 'binary_circulant'
                initial_pattern = randi([0 1], 1,sub_pixels^2);
                cir_temp = gallery('circul', initial_pattern);
                phi = cir_temp(1:m,1:n);

                
         end

    %___THETA___
    %___NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues___
    theta = zeros(m,n);
    for theta_loop = 1:n
        ek = zeros(1,n);
        ek(theta_loop) = 1;
        switch image_transformation_algorithm
            case 'idct'
                psi = idct2(ek)';
            case 'ifwht'
                psi = ifwht(ek)';
        end
        theta(:,theta_loop) = phi*psi;
    end
        frame_index = 1;
        for frame_number = 1:1:10         
            %___LOAD image___
            block_counting = 0;
            frame_number
            load_frame = imread(imageFiles(frame_number).name);
            %load_frame = snapshot(cam);
            if(strcmp(color_mode,'rgb') || strcmp(color_mode,'RGB'))
                frame = load_frame(:,:,:);
                frame = imresize(frame,[720, 1280]);
                plane = 3;
            else
                frame = rgb2gray(load_frame);
                frame = imresize(frame,[720, 1280]);
                plane = 1;
            end
            %___RESET STATE___
            N_1 = zeros(1, size(frame,1)/sub_pixels) + sub_pixels;
            N_2 = zeros(1, size(frame,2)/sub_pixels) + sub_pixels;
            for i = 1:plane
                C(:,:,i) = mat2cell(double(frame(:,:,i)), N_1, N_2);
            end

            %___RESET STATE___
            if(reset == 1)
                reset = 0;
                threshold                   = zeros(size(N_1,2), size(N_2,2), plane);
                quantization.quantization_parameter = zeros(sub_pixels,1);
                y                           = zeros(m,1);
                y_buffer_upleft_encoder     = zeros((m), size(frame,2)/sub_pixels);
                y_buffer_up_encoder         = zeros((m), size(frame,2)/sub_pixels);
                y_buffer_left_encoder       = zeros(m, 1);
                y_buffer_diag_up_encoder    = zeros(m, 1);
                y_buffer_diag_left_encoder  = zeros(m, 1);
                y_buffer_dc_encoder         = zeros(m, 1);
                y_buffer_cp_encoder         = (zeros(m, 1));
                
                y_buffer_up_decoder         = zeros((m), size(frame,2)/sub_pixels);
                y_buffer_upleft_decoder     = zeros((m), size(frame,2)/sub_pixels);
                y_buffer_left_decoder       = zeros(m, 1);
                y_buffer_diag_up_decoder    = zeros(m, 1);
                y_buffer_diag_left_decoder  = zeros(m, 1);
                y_buffer_dc_decoder         = zeros(m, 1);
                y_buffer_cp_decoder         = (zeros(m, 1));
                y_predicted_encoder         = zeros(m, 1);
                plane_and                   = zeros(1,3);
                image_bpp                   = zeros(1,3);
                image_psnr                  = zeros(1,3);
                image_ssim                  = zeros(1,3);
                y_buffer_cp_encoder(1)      = 32640;
                y_buffer_cp_encoder(2:m)  = 16320;
                y_buffer_cp_decoder(1)      = 32640;
                y_buffer_cp_decoder(2:m)  = 16320;
                if(strcmp(color_mode,'rgb') || strcmp(color_mode,'RGB'))
                    bpp_buffer              = zeros(1,3);
                else
                    bpp_buffer              = zeros(1,1);
                end
                temp_padding                = zeros(size(frame,1)+2, size(frame,2)+2, plane);
                res_temp_padding            = zeros(size(frame,1)+2, size(frame,2)+2, plane);
                buffer                      = cell(size(N_1,2), size(N_2,2), plane);
                buffer_1                    = cell(size(N_1,2), size(N_2,2), plane);
                buffer_2                    = cell(size(N_1,2), size(N_2,2), plane);
                estimate_background_encoder = cell(size(N_1,2), size(N_2,2), plane);
                estimate_background_decoder = cell(size(N_1,2), size(N_2,2), plane);
                y_residual                  = cell(size(N_1,2), size(N_2,2), plane);
                combination_frame           = cell(size(N_1,2), size(N_2,2), plane);
                y_quantized                 = cell(size(N_1,2), size(N_2,2), plane);
                y_dequantized               = cell(size(N_1,2), size(N_2,2), plane);
                reconstructed_image         = cell(size(N_1,2), size(N_2,2), plane);
                res_reconstructed_image     = cell(size(N_1,2), size(N_2,2), plane);
                modes                       = zeros(size(frame,1)/sub_pixels, size(frame,2)/sub_pixels, plane);
                bit_shift                   = cell(size(frame,1)/sub_pixels, size(frame,2)/sub_pixels, plane);
                res_video_buffer            = zeros(size(frame,1), size(frame,2), plane, 100);
%                 video_buffer                = zeros(size(frame,1), size(frame,2), plane, 100);
                quartus_output              = zeros(1200, 1, 128);
                for i = 1:size(N_1, 2)
                    for j = 1:size(N_2, 2)
                        for k = 1:plane
                            buffer{i,j,k}                      = zeros(m, 1);
                            buffer_1{i,j,k}                    = zeros(m, 1);
                            buffer_2{i,j,k}                    = zeros(m, 1);
                            estimate_background_encoder{i,j,k} = zeros(m, 1);
                            estimate_background_decoder{i,j,k} = zeros(m, 1);
                            combination_frame{i,j,k}           = zeros(m, 1);
                            y_residual{i,j,k}                  = zeros(m, 1);
                            y_quantized{i,j,k}                 = zeros(m, 1);
                            y_dequantized{i,j,k}               = zeros(m, 1);
                            
                            reconstructed_image{i,j,k}         = zeros(m, 1);
                            res_reconstructed_image{i,j,k}     = zeros(m, 1);
                        end
                    end
                end
            end

            %___THE RANDOM PROJECTION___
            disp('plane');
            for k = 1:plane
                for i = 1:1:size(frame,1)/sub_pixels
                    for j = 1:1:size(frame,2)/sub_pixels
                       one_block_image(:,:,k) = reshape(C{i,j,k}.',1,[])';
                       y = BCS_encoder(one_block_image(:,:,k), phi);
                       if (strcmp(quartus_interface,'on'))
                           %___TAPPING___
%                            for nf = 1:size(y,1)
%                               outputFile = fopen(['Tokyo_y',num2str(nf-1),'.dat'],'at');
%                               fprintf(outputFile, '%c',dec2bin(y(nf)));
%                               fprintf(outputFile, '\n');
%                               fclose(outputFile);
%                           end
                       else
%                             xmin=10;
%                             xmax=100;
%                             n=63;
%                             random=xmin+rand(1,n)*(xmax-xmin);
%                             random = random';
%                             y(2:64) = (y(1)/2)-random;
                            if (j == 1)
                               [y_buffer_left_encoder, y_buffer_up_encoder(:,j), y_buffer_dc_encoder, y_buffer_cp_encoder, y_predicted_encoder, modes(i,j,k)] = coding_method_first_block(y, ...
                                                                                                                                                                  phi, ...
                                                                                                                                                                  i, ...
                                                                                                                                                                  j, ...
                                                                                                                                                                  sub_pixels, ...
                                                                                                                                                                  m, ...
                                                                                                                                                                  n, ...
                                                                                                                                                                  y_buffer_up_encoder(:,j), ...
                                                                                                                                                                  y_buffer_left_encoder, ...
                                                                                                                                                                  y_buffer_dc_encoder, ...
                                                                                                                                                                 y_buffer_cp_encoder, ...
                                                                                                                                                                  'intra_prediction');

                            else
                               [y_buffer_upleft_encoder(:,j-1),y_buffer_left_encoder, y_buffer_up_encoder(:,j), y_buffer_dc_encoder,y_buffer_diag_up_encoder,y_buffer_diag_left_encoder, y_buffer_cp_encoder, y_predicted_encoder, modes(i,j,k)] = coding_method(y, ...
                                                                                                                                                                  phi, ...
                                                                                                                                                                  i, ...
                                                                                                                                                                  j, ...
                                                                                                                                                                  sub_pixels, ...
                                                                                                                                                                  m, ...
                                                                                                                                                                  n, ...
                                                                                                                                                                  y_buffer_upleft_encoder(:,j-1),...
                                                                                                                                                                  y_buffer_up_encoder(:,j), ...
                                                                                                                                                                  y_buffer_left_encoder, ...
                                                                                                                                                                  y_buffer_dc_encoder, ...
                                                                                                                                                                  y_buffer_diag_up_encoder,...
                                                                                                                                                                  y_buffer_diag_left_encoder,...
                                                                                                                                                                 y_buffer_cp_encoder, ...
                                                                                                                                                                  'intra_prediction');

                           end
                            y_residual{i,j,k} = y - y_predicted_encoder;
                            
                    end
                end
                end  
            end
           %___QUANTIZATION___
            bf_search      = 0;
            bits_shift_arr = [];
            quantization_candidates = [];
            bitrate_limit  = 0.3;
            while 1
                bf_search                            =  bf_search + 1;
               quantization_candidates(bf_search)    = round((max(y_residual{i,j,k}) - min(y_residual{i,j,k}))/2^bf_search);
               if(quantization_candidates(bf_search) == 0)
                   break
               end
               disp('------------------------------------------------------------------------------')
               quantization_candidates               = sort(quantization_candidates);
               
               triangle_bit                          = floor(triangle_th(quantization_candidates, size(quantization_candidates,2)));
               if(~isnan(triangle_bit))
                   bits_shift_arr(end+1)             = round((max(y_residual{i,j,k}) - min(y_residual{i,j,k}))/2^triangle_bit);
               end
            end 
             bits_shift_arr                          = sort(bits_shift_arr);
            for step = 1:1:size(bits_shift_arr,2)
                bpp_buffer(k) = 0;
                for k = 1:plane
                    for i = 1:1:size(frame,1)/sub_pixels
                        for j = 1:1:size(frame,2)/sub_pixels
                            bits_shift_arr(step);
                            bits_shift = bits_shift_arr(step);
                            y_quantized{i,j,k}                   = floor(y_residual{i,j,k}./bits_shift_arr(step));
                            bpp_buffer(k)                        = bpp_buffer(k) + measurement_entropy(y_quantized{i,j,k},size(frame,1)*size(frame,2));
                        end
                    end
                end
                bpp_buffer(k);  
                if (bpp_buffer(k) <= bitrate_limit)
                     break;
                end
            end
                                                %%%%%%%%%%%%%%%%%%%%
                                                %%%%% CHANNALS %%%%%
                                                
                                                %%%%%%%%%%%%%%%%%%%%
            if(strcmp(quartus_interface,'on'))
                for file_loop = 1:128
                 quartus_output(:,:,file_loop) = load(strcat('C:/Users/Jay/Desktop/Research/FPGA_ASIC/gray_compressed_sensing_walsh_ifwht_intra_prediction/y_reconstruct/Tokyo_y',num2str(file_loop-1),'.dat'));
                end
                 modes = load(strcat('C:/Users/Jay/Desktop/Research/FPGA_ASIC/gray_compressed_sensing_walsh_ifwht_intra_prediction/y_reconstruct/modes.dat'));
                 modes = reshape(modes.',40,[])';
            end
            for k = 1:plane
                internal_slot=1;
                for i = 1:size(frame,1)/sub_pixels
                    for j = 1:size(frame,2)/sub_pixels
                       if(strcmp(quartus_interface,'on'))
                           y_quantized{i,j,k} = [quartus_output(internal_slot,:,1), quartus_output(internal_slot,:,2), quartus_output(internal_slot,:,3), quartus_output(internal_slot,:,4), quartus_output(internal_slot,:,5), quartus_output(internal_slot,:,6), quartus_output(internal_slot,:,7), quartus_output(internal_slot,:,8), quartus_output(internal_slot,:,9), quartus_output(internal_slot,:,10), ...
                                                 quartus_output(internal_slot,:,11), quartus_output(internal_slot,:,12), quartus_output(internal_slot,:,13), quartus_output(internal_slot,:,14), quartus_output(internal_slot,:,15), quartus_output(internal_slot,:,16), quartus_output(internal_slot,:,17), quartus_output(internal_slot,:,18), quartus_output(internal_slot,:,19), quartus_output(internal_slot,:,20), ...
                                                 quartus_output(internal_slot,:,21), quartus_output(internal_slot,:,22), quartus_output(internal_slot,:,23), quartus_output(internal_slot,:,24), quartus_output(internal_slot,:,25), quartus_output(internal_slot,:,26), quartus_output(internal_slot,:,27), quartus_output(internal_slot,:,28), quartus_output(internal_slot,:,29), quartus_output(internal_slot,:,30), ...
                                                 quartus_output(internal_slot,:,31), quartus_output(internal_slot,:,32), quartus_output(internal_slot,:,33), quartus_output(internal_slot,:,34), quartus_output(internal_slot,:,35), quartus_output(internal_slot,:,36), quartus_output(internal_slot,:,37), quartus_output(internal_slot,:,38), quartus_output(internal_slot,:,39), quartus_output(internal_slot,:,40), ...
                                                 quartus_output(internal_slot,:,41), quartus_output(internal_slot,:,42), quartus_output(internal_slot,:,43), quartus_output(internal_slot,:,44), quartus_output(internal_slot,:,45), quartus_output(internal_slot,:,46), quartus_output(internal_slot,:,47), quartus_output(internal_slot,:,48), quartus_output(internal_slot,:,49), quartus_output(internal_slot,:,50), ...
                                                 quartus_output(internal_slot,:,51), quartus_output(internal_slot,:,52), quartus_output(internal_slot,:,53), quartus_output(internal_slot,:,54), quartus_output(internal_slot,:,55), quartus_output(internal_slot,:,56), quartus_output(internal_slot,:,57), quartus_output(internal_slot,:,58), quartus_output(internal_slot,:,59), quartus_output(internal_slot,:,60), ...
                                                 quartus_output(internal_slot,:,61), quartus_output(internal_slot,:,62), quartus_output(internal_slot,:,63), quartus_output(internal_slot,:,64), quartus_output(internal_slot,:,65), quartus_output(internal_slot,:,66), quartus_output(internal_slot,:,67), quartus_output(internal_slot,:,68), quartus_output(internal_slot,:,69), quartus_output(internal_slot,:,70), ...
                                                 quartus_output(internal_slot,:,71), quartus_output(internal_slot,:,72), quartus_output(internal_slot,:,73), quartus_output(internal_slot,:,74), quartus_output(internal_slot,:,75), quartus_output(internal_slot,:,76), quartus_output(internal_slot,:,77), quartus_output(internal_slot,:,78), quartus_output(internal_slot,:,79), quartus_output(internal_slot,:,80), ...
                                                 quartus_output(internal_slot,:,81), quartus_output(internal_slot,:,82), quartus_output(internal_slot,:,83), quartus_output(internal_slot,:,84), quartus_output(internal_slot,:,85), quartus_output(internal_slot,:,86), quartus_output(internal_slot,:,87), quartus_output(internal_slot,:,88), quartus_output(internal_slot,:,89), quartus_output(internal_slot,:,90), ...
                                                 quartus_output(internal_slot,:,91), quartus_output(internal_slot,:,92), quartus_output(internal_slot,:,93), quartus_output(internal_slot,:,94), quartus_output(internal_slot,:,95), quartus_output(internal_slot,:,96), quartus_output(internal_slot,:,97), quartus_output(internal_slot,:,98), quartus_output(internal_slot,:,99), quartus_output(internal_slot,:,100), ...
                                                 quartus_output(internal_slot,:,101), quartus_output(internal_slot,:,102), quartus_output(internal_slot,:,103), quartus_output(internal_slot,:,104), quartus_output(internal_slot,:,105), quartus_output(internal_slot,:,106), quartus_output(internal_slot,:,107), quartus_output(internal_slot,:,108), quartus_output(internal_slot,:,109), quartus_output(internal_slot,:,110), ...
                                                 quartus_output(internal_slot,:,111), quartus_output(internal_slot,:,112), quartus_output(internal_slot,:,113), quartus_output(internal_slot,:,114), quartus_output(internal_slot,:,115), quartus_output(internal_slot,:,116), quartus_output(internal_slot,:,117), quartus_output(internal_slot,:,118), quartus_output(internal_slot,:,119), quartus_output(internal_slot,:,120), ...
                                                 quartus_output(internal_slot,:,121), quartus_output(internal_slot,:,122), quartus_output(internal_slot,:,123), quartus_output(internal_slot,:,124), quartus_output(internal_slot,:,125), quartus_output(internal_slot,:,126), quartus_output(internal_slot,:,127), quartus_output(internal_slot,:,128)]';
                           for q_i = 1:bits_shift(bits_shift_loop)
                               y_quantized{i,j,k} = bitshift(y_quantized{i,j,k},1);
                           end
                          
                           y_dequantized{i,j,k} = double(typecast(uint16(bin2dec(dec2bin(y_quantized{i,j,k}))),'int16'));
                       else                           %___DE-QUANTIZATION___
                           y_dequantized{i,j,k}  = (y_quantized{i,j,k}*bits_shift) + 0.5;
%                             y_dequantized{i,j,k}                 	   = floor((y_quantized{i,j,k}*quantization.quantization_parameter(i,j,k)) + y_predicted_encoder);  

                       end
                       
                       if(j == 1)
                           if(modes(i,j,k) == 0)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_left_decoder;
                           elseif(modes(i,j,k) == 1)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_up_decoder(:,j);
                           elseif(modes(i,j,k) == 2)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_dc_decoder;
                           elseif(modes(i,j,k) == 4)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_cp_decoder;
                           else
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_cp_decoder;
                           end
                        else
                           if(modes(i,j,k) == 0)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_left_decoder;
                           elseif(modes(i,j,k) == 1)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_up_decoder(:,j);
                           elseif(modes(i,j,k) == 2)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_dc_decoder;
                           elseif(modes(i,j,k) == 3)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_upleft_decoder(:,j-1);
                           elseif(modes(i,j,k) == 5)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_diag_up_decoder;
                           elseif(modes(i,j,k) == 6)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_diag_left_decoder;
                           elseif(modes(i,j,k) == 4)
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_cp_decoder;
                           else
                              combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_cp_decoder;
                           end
                       end
                      if(j == 1)
                          [y_buffer_left_decoder, y_buffer_up_decoder(:,j), y_buffer_dc_decoder, y_buffer_cp_decoder] = coding_method_first_block(combination_frame{i,j,k}, ...
                                                                                                                                      phi, ...
                                                                                                                                      i, ...
                                                                                                                                      j, ...
                                                                                                                                      sub_pixels, ...
                                                                                                                                      m, ...
                                                                                                                                      n, ...
                                                                                                                                      y_buffer_up_decoder(:,j), ...
                                                                                                                                      y_buffer_left_decoder, ...
                                                                                                                                      y_buffer_dc_decoder, ...
                                                                                                                                      y_buffer_cp_decoder, ...
                                                                                                                                      'intra_prediction');

                     else
                          [y_buffer_upleft_decoder(:,j-1),y_buffer_left_decoder, y_buffer_up_decoder(:,j), y_buffer_dc_decoder,y_buffer_diag_up_decoder,y_buffer_diag_left_decoder, y_buffer_cp_decoder] = coding_method(combination_frame{i,j,k}, ...
                                                                                                                                      phi, ...
                                                                                                                                      i, ...
                                                                                                                                      j, ...
                                                                                                                                      sub_pixels, ...
                                                                                                                                      m, ...
                                                                                                                                      n, ...
                                                                                                                                      y_buffer_upleft_decoder(:,j-1),...
                                                                                                                                      y_buffer_up_decoder(:,j), ...
                                                                                                                                      y_buffer_left_decoder, ...
                                                                                                                                      y_buffer_dc_decoder, ...
                                                                                                                                      y_buffer_diag_up_decoder,...
                                                                                                                                      y_buffer_diag_left_decoder,...
                                                                                                                                      y_buffer_cp_decoder, ...
                                                                                                                                      'intra_prediction');
                     end
                      reconstructed_image{i,j,k}         = BCS_reconstruction(combination_frame{i,j,k}, ...
                                                                              theta, ...
                                                                              phi, ...
                                                                              image_transformation_algorithm, ...
                                                                              image_reconstruction_algorithm, ...
                                                                              sub_pixels, ...
                                                                              opts);

                      res_reconstructed_image{i,j,k}     = BCS_reconstruction(y_dequantized{i,j,k}, ...
                                                                              theta, ...
                                                                              phi, ...
                                                                              image_transformation_algorithm, ...
                                                                              image_reconstruction_algorithm, ...
                                                                              sub_pixels, ...
                                                                              opts);
                      internal_slot=internal_slot+1;
                    end
                end
            end
            %___FINAL PROCESS ZONE___%
            for k = 1:plane
                if(plane == 3)
                   temp_padding(:,:,k) = padarray(floor(cell2mat(reconstructed_image(:,:,k))),[1 1],'symmetric','both');
                else
                   temp_padding = padarray(floor(cell2mat(reconstructed_image(:,:))),[1 1],'symmetric','both'); %padding
                end

                %___Overlapped Filtering___ Y X
                for i = 2:size(temp_padding,1)-1
                    for j = 2:size(temp_padding,2)-1
                        video_buffer(i-1,j-1,k,frame_number) = floor((temp_padding(i,j,k)+temp_padding(i,j-1,k)+temp_padding(i,j+1,k))/3);
                    end
                end
            end

            %___FINAL PROCESS ZONE___%
            for k = 1:plane
                if(plane == 3)
                   res_temp_padding(:,:,k) = padarray(floor(cell2mat(res_reconstructed_image(:,:,k))),[1 1],'symmetric','both');
                else
                   res_temp_padding = padarray(floor(cell2mat(res_reconstructed_image(:,:))),[1 1],'symmetric','both'); %padding
                end

                %___Overlapped Filtering___ Y X
                for i = 2:size(res_temp_padding,1)-1
                    for j = 2:size(res_temp_padding,2)-1
                        res_video_buffer(i-1,j-1,k,frame_number) = floor((res_temp_padding(i,j,k)+res_temp_padding(i,j-1,k)+res_temp_padding(i,j+1,k))/3);
                    end
                end
            end
            %___RESET FOR NEW FRAME___
            y_buffer_up_encoder         = zeros((m), size(frame,2)/sub_pixels);
            y_buffer_left_encoder       = zeros(m, 1);
            y_buffer_dc_encoder         = zeros(m, 1);
            y_buffer_diag_up_encoder    = zeros(m,1);
            y_buffer_diag_left_encoder    = zeros(m,1);
            y_buffer_upleft_encoder     = zeros((m), size(frame,2)/sub_pixels);
            y_buffer_cp_encoder         = (zeros(m, 1));

            y_buffer_up_decoder         = zeros((m), size(frame,2)/sub_pixels);
            y_buffer_left_decoder       = zeros(m, 1);
            y_buffer_dc_decoder         = zeros(m, 1);
            y_buffer_diag_up_decoder    = zeros(m,1);
            y_buffer_diag_left_decoder    = zeros(m,1);
            y_buffer_upleft_decoder     = zeros((m), size(frame,2)/sub_pixels);
            y_buffer_cp_decoder         = (zeros(m, 1));
            y_buffer_cp_encoder(1)      = 32640;
            y_buffer_cp_encoder(2:m)  = 16320;
            y_buffer_cp_decoder(1)      = 32640;
            y_buffer_cp_decoder(2:m)  = 16320;
            %modes                       = zeros(size(frame,1)/sub_pixels, size(frame,2)/sub_pixels, plane);
            
            %___QUATITATIVE MATRICES___
            %block_per_frame(frame_number) = block_counting;
            image_bpp(frame_number)  = sum(bpp_buffer) + 1/(sub_pixels * sub_pixels);
            image_psnr(frame_number) = psnr(uint8(video_buffer(:,:,:,frame_number)), frame);
            image_ssim(frame_number) = ssim(uint8(video_buffer(:,:,:,frame_number)), frame);
            bpp_buffer               = zeros(k,1);
        end
        %image = medfilt2((mat2gray(uint8(video_buffer(:,:,:,frame_number)))),'indexed');

       

    save(fullfile(strcat('PSNR_WILD_', num2str(sampling_rate),'_BITS_SHIFT_', '.mat')), 'image_psnr');
    save(fullfile(strcat('SSIM_WILD_', num2str(sampling_rate),'_BITS_SHIFT_', '.mat')), 'image_ssim');
    save(fullfile(strcat('BPP_WILD_' , num2str(sampling_rate),'_BITS_SHIFT_', '.mat')), 'image_bpp');
    video_out = VideoWriter(fullfile(strcat('WILD', ...
                                            '_Frame_Skip_', num2str(sampling_rate), ...
                                            '_', measurement_matrix_construction, ...
                                            '_', image_reconstruction_algorithm, ...
                                            '_', image_transformation_algorithm, ...
                                            '_', color_mode, ...
                                            '_Linear_Filter', ...
                                            '.avi')), 'Uncompressed AVI'); %create the video object
    video_out.FrameRate = 30;                                    
    open(video_out); %open the file for writing
    for loop = 1:frame_number-1
       writeVideo(video_out, uint8(video_buffer(:,:,:,loop)));
    end
    close(video_out);

    video_out = VideoWriter(fullfile(strcat('RESIDUAL_image', ...
                                            '_Frame_Skip_', num2str(sampling_rate), ...
                                            '_', measurement_matrix_construction, ...
                                            '_', image_reconstruction_algorithm, ...
                                            '_', image_transformation_algorithm, ...
                                            '_', color_mode, ...
                                            '_Linear_Filter', ...
                                            '.avi')), 'Uncompressed AVI'); %create the video object
    video_out.FrameRate = 30;                                    
    open(video_out); %open the file for writing
    
    for loop = 1:frame_number-1
       writeVideo(video_out, mat2gray(uint8(res_video_buffer(:,:,:,loop))));
    end
    close(video_out);
end
profile report
profile off