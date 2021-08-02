function [SR] = SRcontrol(max_bit,min_bit, frame_high, frame_width, plane, fps, SR_arr)
    %__Maximum bits comsumption per second
    %start with SR = 0.75
    %SR_arr         = [0.75, 0.5, 0.25, 0.125, 0.065];
    Randbandwidth  = (max_bit - min_bit).*rand(100,1) + min_bit;
    available_bit  = randsample(Randbandwidth,1);
    resolution     = frame_high * frame_width;
    %resolutionWblanking = 
    switch plane
        case 3
            consumption_bit = resolution*3*8*fps;
        case 1
            consumption_bit = resolution*3*8*fps;
    end
    ideaSR    = (available_bit/consumption_bit)
    [val,idx] =  min(abs(ideaSR-SR_arr));
    SR        = SR_arr(idx);
    
    
    
    
    
            
            
    
    