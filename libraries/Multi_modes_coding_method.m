function [y_buffer_left, y_buffer_up, y_buffer_dc,y_buffer_uull,y_buffer_lluu,y_buffer_ulul,y_buffer_lulu, y_buffer_cp, y_predicted, modes] = Multi_modes_coding_method(y, phi, i, j, sub_pixels, m, n, y_buffer_up, y_buffer_left, y_buffer_dc, y_buffer_uull,y_buffer_lluu,y_buffer_ulul,y_buffer_lulu, y_buffer_cp, method)
    switch method
        case 'intra_prediction'
            %___SUM OF ABSOLUTE DIFFERENCE___
            SAD_y_left_mode   = sum(imabsdiff(y,y_buffer_left));
            SAD_y_up_mode     = sum(imabsdiff(y,y_buffer_up));
            SAD_y_cp_mode     = sum(imabsdiff(y,y_buffer_cp));
            SAD_y_dc_mode     = sum(imabsdiff(y,y_buffer_dc));
            SAD_y_uull_mode   = sum(imabsdiff(y,y_buffer_uull));
            SAD_y_lluu_mode   = sum(imabsdiff(y,y_buffer_lluu));
            SAD_y_ulul_mode   = sum(imabsdiff(y,y_buffer_ulul));
            SAD_y_lulu_mode   = sum(imabsdiff(y,y_buffer_lulu));
%  
%             SAD_y_left_mode   = norm(y - y_buffer_left);
%             SAD_y_up_mode     = norm(y - y_buffer_up);
%             SAD_y_cp_mode     = norm(y - y_buffer_cp);
%             SAD_y_dc_mode     = norm(y - y_buffer_dc);
%             SAD_y_uull_mode   = norm(y - y_buffer_uull);
%             SAD_y_lluu_mode   = norm(y - y_buffer_lluu);
%             SAD_y_ulul_mode   = norm(y - y_buffer_ulul);
%             SAD_y_lulu_mode   = norm(y - y_buffer_lulu);
            SAD_candidate     = min([SAD_y_left_mode SAD_y_dc_mode SAD_y_cp_mode SAD_y_up_mode SAD_y_uull_mode SAD_y_lluu_mode SAD_y_ulul_mode SAD_y_lulu_mode]);
            if(SAD_candidate     == SAD_y_dc_mode)
                y_predicted      =  y_buffer_dc;
                modes = 0;
            elseif(SAD_candidate == SAD_y_left_mode)
                y_predicted      =  y_buffer_left;
                modes = 1;
            elseif(SAD_candidate == SAD_y_up_mode)
                y_predicted      =  y_buffer_up;
                modes = 2;
            elseif(SAD_candidate == SAD_y_uull_mode)
                y_predicted      =  y_buffer_uull;
                modes = 3;
            elseif(SAD_candidate == SAD_y_lulu_mode)
                y_predicted      =  y_buffer_lulu;
                modes = 4;
            elseif(SAD_candidate == SAD_y_lluu_mode)
                y_predicted      =  y_buffer_lluu;
                modes = 5;
            elseif(SAD_candidate == SAD_y_ulul_mode)
                y_predicted      =  y_buffer_ulul;
                modes = 6;
            elseif(SAD_candidate == SAD_y_cp_mode)
                y_predicted      =  y_buffer_cp;
                modes = 7;
            else
                y_predicted      =  y_buffer_cp;
                modes = 7;
                
            end
%             if((SAD_y_left_mode < SAD_y_up_mode) && (SAD_y_left_mode < SAD_y_dc_mode) && (SAD_y_left_mode < SAD_y_cp_mode)&& (SAD_y_left_mode < SAD_y_uull_mode)&& (SAD_y_left_mode < SAD_y_lulu_mode) && (SAD_y_left_mode < SAD_y_lluu_mode)&& (SAD_y_left_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_left;
%                 modes = 1;
%             elseif((SAD_y_up_mode < SAD_y_left_mode) && (SAD_y_up_mode < SAD_y_dc_mode) && (SAD_y_up_mode < SAD_y_cp_mode)&& (SAD_y_up_mode < SAD_y_uull_mode)&& (SAD_y_up_mode < SAD_y_lulu_mode)&& (SAD_y_up_mode < SAD_y_lluu_mode)&& (SAD_y_up_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_up;
%                 modes = 2;
%             elseif((SAD_y_dc_mode < SAD_y_up_mode) && (SAD_y_dc_mode < SAD_y_left_mode) && (SAD_y_dc_mode < SAD_y_cp_mode)&& (SAD_y_dc_mode < SAD_y_uull_mode)&& (SAD_y_dc_mode < SAD_y_lulu_mode)&& (SAD_y_dc_mode < SAD_y_lluu_mode)&& (SAD_y_dc_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_dc;
%                 modes = 0;
%             elseif ((SAD_y_uull_mode < SAD_y_up_mode) && (SAD_y_uull_mode < SAD_y_left_mode) && (SAD_y_uull_mode < SAD_y_dc_mode) && (SAD_y_uull_mode < SAD_y_cp_mode)&& (SAD_y_uull_mode < SAD_y_lulu_mode)&& (SAD_y_uull_mode < SAD_y_lluu_mode)&& (SAD_y_uull_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_uull;
%                 modes = 3;
%             elseif ((SAD_y_lulu_mode < SAD_y_up_mode) && (SAD_y_lulu_mode < SAD_y_left_mode) && (SAD_y_lulu_mode < SAD_y_dc_mode) && (SAD_y_lulu_mode < SAD_y_cp_mode)&& (SAD_y_lulu_mode < SAD_y_uull_mode)&& (SAD_y_lulu_mode < SAD_y_lluu_mode)&& (SAD_y_lulu_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_lulu;
%                 modes = 4;
%             elseif ((SAD_y_lluu_mode < SAD_y_up_mode) && (SAD_y_lluu_mode < SAD_y_left_mode) && (SAD_y_lluu_mode < SAD_y_dc_mode) && (SAD_y_lluu_mode < SAD_y_cp_mode)&& (SAD_y_lluu_mode < SAD_y_lulu_mode)&& (SAD_y_lluu_mode < SAD_y_uull_mode)&& (SAD_y_lluu_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_lluu;
%                 modes = 5;
%             elseif ((SAD_y_ulul_mode < SAD_y_up_mode) && (SAD_y_ulul_mode < SAD_y_left_mode) && (SAD_y_ulul_mode < SAD_y_dc_mode) && (SAD_y_ulul_mode < SAD_y_cp_mode)&& (SAD_y_ulul_mode < SAD_y_lulu_mode)&& (SAD_y_ulul_mode < SAD_y_lluu_mode)&& (SAD_y_ulul_mode < SAD_y_uull_mode))
%                 y_predicted   = y_buffer_ulul;
%                 modes = 6;
%             elseif((SAD_y_cp_mode < SAD_y_up_mode) && (SAD_y_cp_mode < SAD_y_left_mode) && (SAD_y_cp_mode < SAD_y_dc_mode)&& (SAD_y_cp_mode < SAD_y_uull_mode)&& (SAD_y_cp_mode < SAD_y_lulu_mode)&& (SAD_y_cp_mode < SAD_y_lluu_mode)&& (SAD_y_cp_mode < SAD_y_ulul_mode))
%                 y_predicted   = y_buffer_cp;
%                 modes = 7;
%             else
%                 y_predicted   = y_buffer_cp;
%                 disp("hi")
%                 modes = 7;
%             end

            coefficient_matrix = sum(phi')';

            %___NEXT ROUND INTRA PREDICTION___
            switch n
                case 256
                        y_buffer_up                = (abs(y(1)-y(2))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_left              = (abs(y(1)-y(32))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_dc                = (y(1)/((sub_pixels*sub_pixels))*coefficient_matrix);
                        y_buffer_lluu              = (y(31)/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_lulu              = (abs(y(1)-y(31))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_uull              = (y(31)/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        y_buffer_ulul              = (abs(y(1)-y(31))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
                        
%                         y_buffer_left_tl_br        = (y(31)/(sub_pixels*sub_pixels/2))*coefficient_matrix;
%                         y_buffer_left_tr_bl        = (abs(y(1)-y(31))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
%                         y_buffer_up_tl_br          = (y(31)/(sub_pixels*sub_pixels/2))*coefficient_matrix;
%                         y_buffer_up_tr_bl          = (abs(y(1)-y(31))/(sub_pixels*sub_pixels/2))*coefficient_matrix;
%                         for i = 1:16:128
%                             if (i <= 64)
%                                 y_buffer_ulul_temp(i:i+7)    = y_buffer_up_tl_br(i:i+7);
%                                 y_buffer_ulul_temp(i+8:i+15) = y_buffer_up_tr_bl(i+8:i+15);     
%                             else
%                                 y_buffer_ulul_temp(i:i+7)    = y_buffer_left_tl_br(i:i+7);
%                                 y_buffer_ulul_temp(i+8:i+15) = y_buffer_left_tl_br(i+8:i+15);
%                             end
%                         end
                        
%                         y_buffer_uull_temp(1:64)   = y_buffer_up(1:64);
%                         y_buffer_uull_temp(65:128) = y_buffer_left(65:128);
%                         y_buffer_uull              = y_buffer_uull_temp.';
%                         y_SB_upper_left                = (abs(y(32)-y(2))/(sub_pixels*sub_pixels/4))*coefficient_matrix(1:32);
%                         y_SB_upper_right               = (abs(y(1)-y(2)-y(32))/(sub_pixels*sub_pixels/4))*coefficient_matrix(33:64);
%                         y_SB_top_left                  = (abs(y(2)-y(32))/(sub_pixels*sub_pixels/4))*coefficient_matrix(33:64);
%                         y_buffer_uull                  = [y_SB_upper_left' y_SB_upper_right' y_SB_top_left' y_SB_upper_right']';
%                         y_buffer_lluu_temp(1:64)   = y_buffer_left(1:64);
%                         y_buffer_lluu_temp(65:128) = y_buffer_up(65:128);
%                         y_buffer_lluu              = y_buffer_lluu_temp.';
%                         for i = 1:16:128
% %                             y_buffer_lulu_temp(i:i+7)    = y_buffer_left(i:i+7);
% %                             y_buffer_lulu_temp(i+8:i+15) = y_buffer_up(i+8:i+15);
%                             y_buffer_ulul_temp(i:i+7)    = y_buffer_up(i:i+7);
%                             y_buffer_ulul_temp(i+8:i+15) = y_buffer_left(i+8:i+15);
%                         end
% %                         y_buffer_lulu                    = y_buffer_lulu_temp.';
%                         y_buffer_ulul                    = y_buffer_ulul_temp.';

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