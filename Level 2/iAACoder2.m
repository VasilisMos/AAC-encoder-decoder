function x = iAACoder2(AACSeq2, fNameOut)

    K = max(size(AACSeq2));   %The number of encoded frames
    x = zeros(K*1024+1024,2); %The output bitstream
    winType = AACSeq2(1).winType; %Assume winType does not change while decoding, as the assignment declares

    s_start = 1-1024;
    s_end = 1024;

%%Extract the time samples for each overlapping frame
    for i=1:K-1
        s_start = s_start+1024;
        s_end = s_end+1024;
        
        frameType = AACSeq2(i).frameType;
        
        switch frameType 
        case "ESH"
            frameF = zeros(128,16);

            frameF(1:128,[1:2:15]) = iTNS(AACSeq2(i).chl.frameF,frameType,AACSeq2(i).chl.TNScoeffs);
            frameF(1:128,[2:2:16]) = iTNS(AACSeq2(i).chr.frameF,frameType,AACSeq2(i).chr.TNScoeffs);

            frameT = iFilterbank(frameF, "ESH", AACSeq2(i).winType);

            x(s_start:s_end,1) = x(s_start:s_end,1) + frameT(:,1);
            x(s_start:s_end,2) = x(s_start:s_end,2) + frameT(:,2);
        otherwise
            frameFchl = iTNS(AACSeq2(i).chl.frameF,frameType,AACSeq2(i).chl.TNScoeffs);
            frameFchr = iTNS(AACSeq2(i).chr.frameF,frameType,AACSeq2(i).chr.TNScoeffs);

            frameF = [frameFchl frameFchr];
            frameT = iFilterbank(frameF, frameType, winType);

            x(s_start:s_end,1) = x(s_start:s_end,1) + frameT(:,1);
            x(s_start:s_end,2) = x(s_start:s_end,2) + frameT(:,2);
        end
    end

    %Write the decoded waveform on the output file (48kHz sampling rate)
    audiowrite(fNameOut,x,48000);

end