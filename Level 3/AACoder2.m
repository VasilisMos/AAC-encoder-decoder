function AACSeq2 = AACoder2(fNameIn)
    
    x = padinput(audioread(fNameIn),1);
    N = max(size(x)); winType = "SIN";

%%Break the samples into overlapping frames of width 2048 each
    frames = []; K = 0;
    inds = 1-1024:1024;

    while not( inds(end) == N )
        inds = inds + 1024;
        K = K+1;

        frames = [frames x(inds,:);];
    end
    
%%Sequence segmentation control 
    AACSeq2(1).frameType = "OLS";  AACSeq2(1).winType = winType;
    AACSeq2(K).frameType = "OLS";  AACSeq2(K).winType = winType;

    for i=2:K-1
        inds = [i*2-1 i*2];
        currFr = frames(:,inds);
        nextFr = frames(:,inds+2);
        AACSeq2(i).frameType = SSC( currFr, nextFr ,AACSeq2(i-1).frameType);
        AACSeq2(i).winType = winType;
    end

%%MDCT for each frame

    for i=1:K-1
        %disp(i)
        frame = frames(:,[i*2-1 i*2]);
        frameF = filterbank(frame, AACSeq2(i).frameType, AACSeq2(i).winType);

        switch AACSeq2(i).frameType
        case "ESH"
            [AACSeq2(i).chl.frameF, AACSeq2(i).chl.TNScoeffs] = TNS(reshape(frameF(:,1),[128 8]),AACSeq2(i).frameType);
            [AACSeq2(i).chr.frameF, AACSeq2(i).chr.TNScoeffs] = TNS(reshape(frameF(:,2),[128 8]),AACSeq2(i).frameType);
        otherwise
            [AACSeq2(i).chl.frameF, AACSeq2(i).chl.TNScoeffs] = TNS(frameF(:,1),AACSeq2(i).frameType);
            [AACSeq2(i).chr.frameF, AACSeq2(i).chr.TNScoeffs] = TNS(frameF(:,2),AACSeq2(i).frameType);
        end
    end
end