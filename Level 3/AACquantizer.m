function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
    
    tbl = load('TableB219.mat');
    MQ =8191;
    
    switch frameType
    case "ESH"
        tbl = tbl.B219b;
        subframes = 8;
        Nb = max(size(tbl));
        N = 128;
    otherwise
        tbl = tbl.B219a;
        subframes = 1;
        Nb = max(size(tbl));
        N=1024;
    end

    S = zeros(N,subframes);
    sfc = zeros(Nb,subframes);

    w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
    w_high = tbl(:,3) + 1; %+1 since matlab is 1-based

    %Calculate 1st scalefactor approximation
    for i=1:subframes
        a_aprox(i) = 16/3*log2(max(frameF(:,i))^(3/4)/MQ);
    end
    
    %Calculate T threshold
    P = zeros(Nb, subframes);
    Perror = zeros(Nb, subframes);
    for b = 1:Nb
        wl = w_low(b);  wh = w_high(b);
        P(b, :) = sum(frameF(wl:wh, :).^2, 1);
    end

    T = P./SMR;
    a = a_aprox(1)*ones(Nb,subframes);
    temp_a = a;
    
    
    %Iterative proccess to find scalefactors and compute S
    for i = 1:subframes
        flags = ones(Nb,1);

        while (max(abs(diff(temp_a(:,i)))) <= 60) 
            a(:,i) = temp_a(:,i);
            
            for b = 1:Nb
                if flags(b)
                    wl = w_low(b);
                    wh = w_high(b);

                    S(wl:wh, i) = quantV3(frameF(wl:wh,i),a(b,i));

                    %Quantization Error
                    Perror(b,i) = sum((frameF(wl:wh,i) - iquantV3(S(wl:wh,i),a(b,i))).^2);
                end
            end
            %Checks if error is above the T threshold
            flags = (T(:, i) > Perror(:,i));
            temp_a(:,i) = a(:,i);

            %for those bands that are not above T threshold, increase its scalefactor unitarily
            temp_a(flags,i) = temp_a(flags,i) + 1;

            if not(any(flags)) %All b's have T(b) < Quant Error
                break;
            end
        end
    end

    %Transform S in the correct format
    S = reshape(S,[1024 1]);
    
    %DPCM To get the scalefactor values
    sfc = a(1, :) - a;
    sfc(2:end, :) = sfc(2:end, :) - sfc(1:end-1, :); 

    %Keep the a(1,:) as Gains
    G = a(1, :);
    
end

function q = quantV3(F,a)
    
    magicNumber = 0.4054;
    MQ = 8191;

    q = sign(F) .* floor((abs(F) .* 2^(-0.25*a) ).^(3/4)  + magicNumber);

end

function x = iquantV3(q,a)
    x = sign(q) .* (abs(q).^(4/3)) .* 2.^(1/4 * a);
end