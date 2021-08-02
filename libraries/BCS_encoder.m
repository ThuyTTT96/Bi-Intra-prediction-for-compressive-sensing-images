function y = BCS_encoder(one_block_image, measurement_matrix)
    %___COMPRESSION___
    y = measurement_matrix*one_block_image;
end
