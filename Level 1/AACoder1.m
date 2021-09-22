function AACSeq1 = AACoder1(fNameIn)
    
    x = padinput(audioread(fNameIn));
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
    AACSeq1(1).frameType = "OLS";  AACSeq1(1).winType = winType;
    AACSeq1(K).frameType = "OLS";  AACSeq1(K).winType = winType;

    for i=2:K-1
        inds = [i*2-1 i*2];
        currFr = frames(:,inds);
        nextFr = frames(:,inds+2);
        AACSeq1(i).frameType = SSC( currFr, nextFr ,AACSeq1(i-1).frameType);
        AACSeq1(i).winType = winType;
    end

%%MDCT for each frame

    for i=1:K
        frame = frames(:,[i*2-1 i*2]);
        frameF = filterbank(frame, AACSeq1(i).frameType, AACSeq1(i).winType);

        switch AACSeq1(i).frameType
        case "ESH"
            AACSeq1(i).chl.frameF = reshape(frameF(:,1),[128 8]);
            AACSeq1(i).chr.frameF = reshape(frameF(:,2),[128 8]);
        otherwise
            AACSeq1(i).chl.frameF = frameF(:,1);
            AACSeq1(i).chr.frameF = frameF(:,2);
        end
    end
end