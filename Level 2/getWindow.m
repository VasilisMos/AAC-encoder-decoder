function win = getWindow(frameType,winType)

    switch winType 
    case "KBD"
        switch frameType
        case "OLS"
            alpha = 6;
            win = getKaiser(2048,alpha);
        case "LSS"
            win = lssKBD();
        case "LPS"
            win = lpsKBD();
        case "ESH"
            alpha = 4;
            win = getKaiser(256,alpha);
        end
    case "SIN"
        switch frameType
        case "OLS"
            N = 2048;
            win = sin(pi/N * ((1:N) - 0.5))';
        case "LSS"
            win = lssSIN();
        case "LPS"
            win = lpsSIN();
        case "ESH"
            N = 256;
            win = sin(pi/N * ((1:N) - 0.5))';
        end
    otherwise
        fprintf("INVALID WINDOW NAME (getWindow())\n");
    end

end

function win = getKaiser(N,alpha)

    kz = kaiser(N/2+1,alpha);
    scale = sum(kz); 
    win = zeros(N,1);

    for i=1:N/2
        win(i) = sum(kz(1:i));
    end

    %the window is symmetric 
    win(N/2+1:end) = win(N/2:-1:1);

    win = sqrt(win/scale);
end

function win = lssKBD()
    N=2048; longA = 6; shortA = 4;
    temp = getKaiser(N,longA);     temp = temp(1:N/2);
    temp2 = getKaiser(256,shortA);  temp2 = temp2(129:end);
    winR = [temp' ones(1,448) temp2'  zeros(1,448)];

    win = winR';
end

function win = lpsKBD()
    win = lssKBD();
    win = win(end:-1:1);
end

function win = lssSIN()
    N = 2048;
    temp = sin(pi/N * ((1:N) - 0.5));    temp = temp(1:N/2);
    temp2 = sin(pi/256 * ((1:256) - 0.5));  temp2 = temp2(129:end);
    winR = [temp ones(1,448) temp2  zeros(1,448)];

    win = winR';

end

function win = lpsSIN()
    win = lssSIN();
    win = win(end:-1:1);
end