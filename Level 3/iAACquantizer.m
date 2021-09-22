function frameF = iAACquantizer(S, sfc, G, frameType)
    
    tbl = load('TableB219.mat');
    
    switch frameType 
    case "ESH"
        tbl = tbl.B219b;
        w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
        w_high = tbl(:,3) + 1; %+1 since matlab is 1-based
        
        Nb = max(size(tbl));
        subframes = 8;
        N = 128;

        S = reshape(S, [N subframes]);
        sfc = reshape(sfc, [42 8]);

        frameF = zeros(N, subframes);
    otherwise
        tbl = tbl.B219a;
        w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
        w_high = tbl(:,3) + 1; %+1 since matlab is 1-based
        
        Nb = max(size(tbl));
        subframes = 1;
        N = 1024;

        frameF = zeros(N, subframes);
    end  

    for ch = 1:subframes
        %Extract scalefactors from DPCM
        a = zeros(Nb,1); a(1) = G(ch);
        for j=2:Nb
            a(j) =  a(1) - sum(sfc(1:j,ch));
        end
        
        %Inverse quantize
        for b = 1:Nb
            wl = w_low(b);
            wh = w_high(b);
            
            frameF(wl:wh, ch) = sign(S(wl:wh, ch)) .* (abs(S(wl:wh, ch)).^(4/3)) .* 2.^(0.25*a(b));              
        end
    end
end