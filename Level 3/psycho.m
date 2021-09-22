function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)
    N = max(size(frameT));
    tbl = load('TableB219.mat');
    bvalLong = tbl.B219a(:,5);
    bvalShort = tbl.B219b(:,5);
    spreads = initspreads(tbl,frameType);

    NMT = 6;
    TMN = 18;
    TBCLIP = 1;

    switch frameType
    case "ESH"
        Nb = length(bvalShort); N = 256;
        tbl = tbl.B219b;
        SMR = zeros(Nb,8);

        %Hann Window
        win = 0.5-0.5*cos(pi/N*((0:N-1)+0.5))';
        
        ind_beg = 449;
        ind_end = 704;
        for ch=1:8
            sw = frameT(ind_beg:ind_end).*win;

            if ch == 1
                prev1T = frameTprev1(1345:1600);
                prev2T = frameTprev1(1217:1472);
            elseif ch == 2
                prev1T = frameT(ind_beg-128:ind_end-128);
                prev2T = frameTprev1(1345:1600);
            else
                prev1T = frameT(ind_beg-128:ind_end-128);
                prev2T = frameT(ind_beg-256:ind_end-256);
            end

            ind_beg = ind_beg + 128;
            ind_end = ind_end + 128;
            
            %Calculate fft
            temp_ind = 1:128;

            Sw = fft(sw); rw = abs(Sw(temp_ind)); fw = angle(Sw(temp_ind));
            Swm1 = fft(prev1T.*win); rwm1 = abs(Swm1(temp_ind)); fwm1 = angle(Swm1(temp_ind));
            Swm2 = fft(prev2T.*win); rwm2 = abs(Swm2(temp_ind)); fwm2 = angle(Swm2(temp_ind));

            rpred = 2*rwm1 - rwm2; %128x1
            fpred = 2*fwm1 - fwm2;

            cw =  sqrt(   (rw.*cos(fw) - rpred.*cos(fpred)).^2   +   (rw.*sin(fw) - rpred.*sin(fpred)).^2  ) ./  (rw + abs(rpred));

            w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
            w_high = tbl(:,3) + 1; %+1 since matlab is 1-based

            for j=1:Nb
                wl = w_low(j);
                wh = w_high(j);

                e(j) = sum(rw(wl:wh).^2);
                cb(j) = sum( cw(w_low(j):w_high(j)).*(rw(w_low(j):w_high(j)).^2)  );
            end
            
            for b=1:Nb
                ecb(b) = e*spreads(:,b);
                ct(b) = cb*spreads(:,b);
            end

            cb = ct./ecb;
            en = ecb./(sum(spreads,1));

            tb = -0.299 - 0.43*log(cb); tb';
            if TBCLIP 
                tb = max(tb,0); tb = min(tb,1);
            end

            SNR = tb*TMN + (1-tb)*NMT;
        
            bc = 10.^(-SNR/10);

            nb = en.*bc;

            for b=1:Nb
                qthr = eps()*N/2*10^(tbl(b,6)/10);
                npart(b) = max(nb(b),qthr);
            end
            temp = (e./npart);
            SMR(:,ch) = temp';
        end
    otherwise %Long window
        Nb = length(bvalLong);
        tbl = tbl.B219a;
        
        %Define Hann Window
        win = 0.5-0.5*cos(pi/N*((0:N-1)+0.5))';

        %Calculate fft
        temp_ind = 1:1024;

        Sw = fft(frameT.*win); rw = abs(Sw(temp_ind)); fw = angle(Sw(temp_ind));
        Swm1 = fft(frameTprev1.*win); rwm1 = abs(Swm1(temp_ind)); fwm1 = angle(Swm1(temp_ind));
        Swm2 = fft(frameTprev2.*win); rwm2 = abs(Swm2(temp_ind)); fwm2 = angle(Swm2(temp_ind));

        rpred = 2*rwm1 - rwm2; %1024x1
        fpred = 2*fwm1 - fwm2;

        cw =  sqrt(   (rw.*cos(fw) - rpred.*cos(fpred)).^2   +   (rw.*sin(fw) - rpred.*sin(fpred)).^2  ) ./  (rw + abs(rpred));
        
        w_low =  tbl(:,2) + 1; %+1 since matlab is 1-based
        w_high = tbl(:,3) + 1; %+1 since matlab is 1-based

        for j=1:Nb
            e(j) = sum(rw(w_low(j):w_high(j)).^2);
            cb(j) = sum( cw(w_low(j):w_high(j)).*(rw(w_low(j):w_high(j)).^2)  );
        end

        for b=1:Nb
            ecb(b) = e*spreads(:,b);
            ct(b) = cb*spreads(:,b);
        end

        cb = ct./ecb;
        en = ecb./(sum(spreads,1));

        %To enable tb clipping to (0,1) comment out line 123
        tb = -0.299 - 0.43*log(cb);
        if TBCLIP 
            tb = max(tb,0); tb = min(tb,1);
        end

        SNR = tb*TMN + (1-tb)*NMT;
        
        bc = 10.^(-SNR/10);

        nb = en.*bc;

        for b=1:Nb
            qthr = eps()*N/2*10^(tbl(b,6)/10);
            npart(b) = max(nb(b),qthr);
        end

        SMR = (e./npart)';
    end

end

function spreads = initspreads(tbl,frameType)

    switch frameType 
    case "ESH"
        bval = tbl.B219b(:,5);
        Nb = length(bval);
        spreads = zeros(Nb,Nb);
    otherwise
        bval = tbl.B219a(:,5);
        Nb = length(bval);
        spreads = zeros(Nb,Nb);
    end

    for i=1:Nb
        for j=1:Nb
            spreads(i,j) = spreadingFun(i,j,bval(i),bval(j));
        end
    end

end

function x = spreadingFun(i,j,bvali,bvalj)
%The spreading Function given as pseudocode in the assignment
    if i>=j
        tmpx = 3.0 * (bvalj - bvali) ;
    else
        tmpx = 1.5 * (bvalj - bvali) ;
    end
    tmpz = 8 * min( (tmpx-0.5)^2 -2*(tmpx-0.5),0);
    tmpy = 15.811389 + 7.5*(tmpx+0.474) - 17.5*sqrt((1+(tmpx+0.474)^2));

    if tmpy<-100
        x = 0;
    else
        x = 10^((tmpz+tmpy)/10);
    end
end


