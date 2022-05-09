function [freq, S_fft] = calc_fft(S, Fs)
    
    for i = 1:size(S,1)
        data = S(i,:);
        L = length(data);
        f = Fs*(0:ceil(L/2))/L;
        Y = fft(data);
        P2 = abs(Y/L);
        P1 = P2(1:ceil(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
    
        
        S_fft(i,:) = P1;
        freq(i,:) = f;
    end  
end

