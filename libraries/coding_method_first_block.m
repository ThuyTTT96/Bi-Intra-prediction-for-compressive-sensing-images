 function [y_buffer_left, y_buffer_up, y_buffer_dc, y_buffer_cp, y_predicted, modes] = coding_method_first_block(y, phi, i, j, sub_pixels, m, n, y_buffer_up, y_buffer_left, y_buffer_dc, y_buffer_cp, method)
    switch method
        case 'intra_prediction'
            % Position of current block
            
            %___SUM OF ABSOLUTE DIFFERENCE___
            SAD_y_left_mode   = sum(imabsdiff(y,y_buffer_left));
            SAD_y_up_mode     = sum(imabsdiff(y,y_buffer_up));
            SAD_y_cp_mode     = sum(imabsdiff(y,y_buffer_cp));
            SAD_y_dc_mode     = sum(imabsdiff(y,y_buffer_dc));
            SAD_candidate     = min([SAD_y_left_mode SAD_y_dc_mode SAD_y_cp_mode SAD_y_up_mode]);

            if (SAD_candidate == SAD_y_left_mode)
                y_predicted   = y_buffer_left;
                modes = 0;
            elseif (SAD_candidate == SAD_y_up_mode)
                y_predicted   = y_buffer_up;
                modes = 1;
            elseif (SAD_candidate == SAD_y_dc_mode)
                y_predicted   = y_buffer_dc;
                modes = 2;
            elseif (SAD_candidate == SAD_y_cp_mode)
                y_predicted   = y_buffer_cp;
                modes = 6;
            else
                y_predicted   = y_buffer_cp;
                modes = 6;
            end

                
            coefficient_matrix = sum(phi')';

            %___NEXT ROUND INTRA PREDICTION___
            switch n
                case 16

                            y_buffer_up           = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_left         = (abs(y(1)-y(8))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_dc           = (y_buffer_up + y_buffer_left)/2;
                case 64

                            y_buffer_up           = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_left         = (abs(y(1)-y(8))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                            y_buffer_dc           = (y_buffer_up + y_buffer_left)/2;                            
                case 256
                        y_buffer_up           = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_left         = (abs(y(1)-y(32))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_dc          = (y_buffer_up + y_buffer_left)/2;

            end

    end
end