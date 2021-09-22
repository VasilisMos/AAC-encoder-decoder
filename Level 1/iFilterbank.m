function frameT = iFilterbank(frameF, frameType, winType)

    switch frameType
    case "ESH"
        alpha = 4;

        frameT = zeros(2048,2);
        indsT = 449:704;

        for i=1:8
            frameT(indsT,1) = frameT(indsT,1) + myIMDCT( frameF(:,i*2-1), winType,alpha,frameType );
            frameT(indsT,2) = frameT(indsT,2) + myIMDCT( frameF(:,i*2), winType,alpha,frameType );
            
            indsT = indsT + 128;
        end
    otherwise
        alpha = 6;
        
        frameT(:,1) = myIMDCT(frameF(:,1),winType,alpha,frameType);
        frameT(:,2) = myIMDCT(frameF(:,2),winType,alpha,frameType);
    end
end


function frameT = myIMDCT(frameF, winType,alpha,frameType)

    frameT = imdct4(frameF)';
    win = getWindow(frameType,winType);

    frameT = win.*(frameT');
end