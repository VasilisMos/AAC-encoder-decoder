%frameFin: MDCT values before TNS (128x8 for short frames, 1024x1 for long frames)
%frameFout: MDCT values after TNS (128x8 for short frames, 1024x1 for long frames)
%TNScoeffs: quantized TNS coefficients (4x8 for ESH, 4x1 else)

function [frameFout, TNScoeffs] = TNS(frameFin, frameType)
    %(IMPORTANT THE FILE MUST BE ON THE SAME FOLDER WITH THE FUNCTIONS)
    %Reads tables from 'TableB219.mat' 
    tbl = load('TableB219.mat');

    

    switch frameType
    case "ESH"  
        tbl = tbl.B219b;

        %Calculate the P(j)
        Nb = size(tbl,1);

        w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
        w_high = tbl(:,3) + 1; %+1 since matlab is 1-based

        for i=1:8
            for j = 1:Nb
                wl = w_low(j);
                wh = w_high(j);
                P(j) = sum(frameFin(wl:wh,i).^2);
                Sw(wl:wh) = sqrt(P(j));
            end
    
            %Smoothing that is presented on the assignment pdf (page 9)
            for k = 127 : -1 : 1
                Sw(k) = (Sw(k) + Sw(k + 1))/2;
            end
            for k = 2 : +1 : 128 
                Sw(k) = (Sw(k) + Sw(k - 1))/2;
            end
    
            Xw = (frameFin(:,i))./(Sw');

            r = autocorr(Xw,4);
            R = toeplitz(r(1:4));
    
            TNScoeffs(:,i) = R\r(2:5);
            TNScoeffs(:,i) = quantize(TNScoeffs(:,i));
            
            frameFout(:,i) = filter([1 -(TNScoeffs(:,i)')],[1],frameFin(:,i));
        end         
    otherwise
        tbl = tbl.B219a;
        Nb = size(tbl,1);

        w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
        w_high = tbl(:,3) + 1; %+1 since matlab is 1-based

        %Calculate P(j) Sw(k)

        for j = 1:Nb
            wl = w_low(j);
            wh = w_high(j); 
            P(j) = sum(frameFin(wl:wh).^2);
            Sw(wl:wh) = sqrt(P(j));
        end

        %Smoothing that is presented on the assignment pdf (page 9)
        for k = 1023 : -1 : 1
            Sw(k) = (Sw(k) + Sw(k + 1))/2;
        end
        for k = 2 : +1 : 1024 
            Sw(k) = (Sw(k) + Sw(k - 1))/2;
        end

        Xw = frameFin./(Sw');

        r = autocorr(Xw,4);
        R = toeplitz(r(1:4));

        TNScoeffs = R\r(2:5);
        TNScoeffs = quantize(TNScoeffs);
        
        frameFout = filter([1 -(TNScoeffs')],[1],frameFin);
    end


end

function q = quantize(x)
    N = length(x);

    levels = -0.75:0.1:0.70;

    for i = 1:N
        if x(i) < levels(1) 
            q(i) = levels(1); continue;
        end
        %disp(i)
        r = (x(i)-levels) >= 0;
        temp = levels(r);
        q(i) = temp(end)+0.05;
        
    end

    q = q';

end