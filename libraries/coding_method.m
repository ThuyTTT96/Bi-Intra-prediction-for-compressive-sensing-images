function [y_buffer_upleft, y_buffer_left, y_buffer_up, y_buffer_dc, y_buffer_diag_up,y_buffer_diag_left, y_buffer_cp, y_predicted, modes] = coding_method(y, phi, i, j, sub_pixels, m, n, y_buffer_upleft,y_buffer_up, y_buffer_left, y_buffer_dc,y_buffer_diag_up, y_buffer_diag_left, y_buffer_cp, method)
    switch method
        case 'intra_prediction'
            % Position of current block
            
            %___SUM OF ABSOLUTE DIFFERENCE___
            SAD_y_left_mode   = sum(imabsdiff(y,y_buffer_left));
            SAD_y_up_mode     = sum(imabsdiff(y,y_buffer_up));
            SAD_y_cp_mode     = sum(imabsdiff(y,y_buffer_cp));
            SAD_y_dc_mode     = sum(imabsdiff(y,y_buffer_dc));
            SAD_y_upleft_mode = sum(imabsdiff(y,y_buffer_upleft));
            SAD_y_diag_up_mode= sum(imabsdiff(y,y_buffer_diag_up));
            SAD_y_diag_left_mode = sum(imabsdiff(y,y_buffer_diag_left));
            SAD_candidate     = min([SAD_y_left_mode SAD_y_dc_mode SAD_y_cp_mode SAD_y_up_mode SAD_y_upleft_mode SAD_y_diag_up_mode SAD_y_diag_left_mode]);

            if (SAD_candidate == SAD_y_left_mode)
                y_predicted   = y_buffer_left;
                modes = 0;
            elseif (SAD_candidate == SAD_y_up_mode)
                y_predicted   = y_buffer_up;
                modes = 1;
            elseif (SAD_candidate == SAD_y_dc_mode)
                y_predicted   = y_buffer_dc;
                modes = 2;
            elseif (SAD_candidate == SAD_y_upleft_mode)
                y_predicted   = y_buffer_upleft;
                modes = 3;
            elseif (SAD_candidate == SAD_y_diag_up_mode)
                y_predicted   = y_buffer_diag_up;
                modes = 5;
            elseif (SAD_candidate == SAD_y_diag_left_mode)
                y_predicted   = y_buffer_diag_left;
                modes = 6;
            elseif (SAD_candidate == SAD_y_cp_mode)
                y_predicted   = y_buffer_cp;
                modes = 4;
            else
                y_predicted   = y_buffer_cp;
                modes = 4;
            end
            coefficient_matrix = sum(phi')';
            %___NEXT ROUND INTRA PREDICTION___
            switch n
                case 16
                            y_buffer_up           = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_left         = (abs(y(1)-y(8))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_upleft       = (abs(2*y(1)-y(2)-y(8))/(sub_pixels*sub_pixels))*coefficient_matrix;
                            y_buffer_dc          = (y_buffer_up + y_buffer_left)/2;
                            y_buffer_diag_up     = (y_buffer_up + y_buffer_upleft)/2;
                            y_buffer_diag_left   = (y_buffer_left + y_buffer_upleft)/2;
                case 256
                            y_buffer_up           = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_left         = (abs(y(1)-y(32))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_upleft       = (abs(2*y(1)-y(2)-y(32))/(sub_pixels*sub_pixels))*coefficient_matrix;
                            y_buffer_dc          = (y_buffer_up + y_buffer_left)/2;
                            y_buffer_diag_up     = (y_buffer_up + y_buffer_upleft)/2;
                            y_buffer_diag_left   = (y_buffer_left + y_buffer_upleft)/2;

                case 128
                            y_buffer_up           = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_left         = (abs(y(1)-y(32))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_upleft       = (abs(2*y(1)-y(2)-y(32))/(sub_pixels*sub_pixels))*coefficient_matrix;
                            y_buffer_dc          = (y_buffer_up + y_buffer_left)/2;
                            y_buffer_diag_up     = (y_buffer_up + y_buffer_upleft)/2;
                            y_buffer_diag_left   = (y_buffer_left + y_buffer_upleft)/2;

            end
            
        case 'motion_detect'
           foreground_searching   = floor(immse(y, estimate_background_encoder{i,j,k}));
           %___MEASUREMENT DECOMPOSITION___
           if(foreground_searching >= threshold(i,j,k))
               foreground{i,j,k} = y;
           else
               foreground{i,j,k} = zeros(m, 1);
           end    
           buffer_2{i,j,k}   = buffer_1{i,j,k};
           buffer_1{i,j,k}   = buffer{i,j,k};
           temp_q            = floor(bitsra(y, bits_shift(bits_shift_loop)));
           temp_deq          = floor(bitsll(temp_q, bits_shift(bits_shift_loop)));
           buffer{i,j,k}     = temp_deq;
           x_1               = immse(estimate_background_encoder{i,i,k},buffer{i,j,k});
           x_2               = immse(estimate_background_encoder{i,j,k},buffer_1{i,j,k});
           x_3               = immse(estimate_background_encoder{i,j,k},buffer_2{i,j,k});
           min_x = min([x_1, x_2, x_3]);
           if(min_x == x_1)
                estimate_background_encoder{i,j,k} = (estimate_background_encoder{i,j,k} + buffer{i,j,k})/2;
           elseif(min_x == x_2)
                estimate_background_encoder{i,j,k} = (estimate_background_encoder{i,j,k} + buffer_1{i,j,k})/2;
           elseif(min_x == x_3)
                estimate_background_encoder{i,j,k} = (estimate_background_encoder{i,j,k} + buffer_2{i,j,k})/2;
           end
           threshold(i,j,k) = median(foreground{i,j,k});
    end
end