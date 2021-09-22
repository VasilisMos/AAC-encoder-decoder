function x = iAACoder1(AACSeq1, fNameOut)

    K = max(size(AACSeq1));   %The number of encoded frames
    x = zeros(K*1024+1024,2); %The output bitstream
    winType = AACSeq1(1).winType; %Assume winType does not change while decoding, as the assignment declares

    s_start = 1-1024;
    s_end = 1024;

%%Extract the time samples for each overlapping frame
    for i=1:K-1
        s_start = s_start+1024;
        s_end = s_end+1024;
        
        frameType = AACSeq1(i).frameType;
        
        switch frameType 
        case "ESH"
            frameF = zeros(128,16);

            frameF(1:128,[1:2:15]) = AACSeq1(i).chl.frameF;
            frameF(1:128,[2:2:16]) = AACSeq1(i).chr.frameF;

            frameT = iFilterbank(frameF, "ESH", AACSeq1(i).winType);

            x(s_start:s_end,1) = x(s_start:s_end,1) + frameT(:,1);
            x(s_start:s_end,2) = x(s_start:s_end,2) + frameT(:,2);
        otherwise
            frameF = [AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF];
            frameT = iFilterbank(frameF, frameType, winType);

            x(s_start:s_end,1) = x(s_start:s_end,1) + frameT(:,1);
            x(s_start:s_end,2) = x(s_start:s_end,2) + frameT(:,2);
        end
    end

    %Write the decoded waveform on the output file (48kHz sampling rate)
    audiowrite(fNameOut,x,48000);

end