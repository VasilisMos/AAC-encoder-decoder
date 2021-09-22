function frameF = filterbank(frameT, frameType, winType)

    switch frameType
    case "ESH"
    
        alpha = 4;

        frameF = zeros(128,16);
        indsT = 449:704;

        for i =1:8
            cols = [(i-1)*2+1 (i-1)*2+2];
            frameF(:,cols(1)) = mymdct(frameT(indsT,1),alpha,winType,frameType);
            frameF(:,cols(2)) = mymdct(frameT(indsT,2),alpha,winType,frameType);

            indsT = indsT + 128;
        end

        leftF = frameF(:,1:2:15);
        rightF = frameF(:,2:2:16); 

        %Make the frameF 1024x2
        frameF = [leftF(:) rightF(:)];

    otherwise 
        alpha = 6;

        frameF(:,1) = mymdct(frameT(:,1),alpha,winType,frameType);
        frameF(:,2) = mymdct(frameT(:,2),alpha,winType,frameType);
    
    end
end


function X = mymdct(x,alpha,winType,frameType)
    N = length(x);
    X = zeros(N/2,1);

    win = getWindow(frameType,winType);

    xW = x.*win;
    X = mdct4(xW);

end