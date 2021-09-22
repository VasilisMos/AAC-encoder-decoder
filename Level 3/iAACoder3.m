function x = iAACoder3(AACSeq3, fNameOut)

    K = max(size(AACSeq3));   %The number of encoded frames
    x = zeros(K*1024+1024,2); %The output bitstream
    winType = AACSeq3(1).winType; %Assume winType does not change while decoding, as the assignment declares

    huffLUT = loadLUT();

    s_start = 1-1024;
    s_end = 1024;

%%Extract the time samples for each overlapping frame
    for i=1:K-1
        s_start = s_start+1024;
        s_end = s_end+1024;
        
        frameType = AACSeq3(i).frameType;
        
        switch frameType 
        case "ESH"
            frameF = zeros(128,16);
            sfcl = zeros(42,8);
            sfcr = zeros(42,8);

            sfcl = decodeHuff(AACSeq3(i).chl.sfc, 12, huffLUT);
            sfcr = decodeHuff(AACSeq3(i).chr.sfc, 12, huffLUT);

            sfcl = reshape(sfcl,[42 8]);
            sfcr = reshape(sfcr,[42 8]);
            
            Sl = decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook, huffLUT);
            Sr = decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook, huffLUT);

            frameFL = iAACquantizer(Sl, sfcl, AACSeq3(i).chl.G, frameType);
            frameFR = iAACquantizer(Sr, sfcr, AACSeq3(i).chr.G, frameType);

            frameF(1:128,[1:2:15]) = iTNS(frameFL,frameType,AACSeq3(i).chl.TNScoeffs);
            frameF(1:128,[2:2:16]) = iTNS(frameFR,frameType,AACSeq3(i).chr.TNScoeffs);

            frameT = iFilterbank(frameF, "ESH", AACSeq3(i).winType);

            x(s_start:s_end,1) = x(s_start:s_end,1) + frameT(:,1);
            x(s_start:s_end,2) = x(s_start:s_end,2) + frameT(:,2);
        otherwise

            sfcl = decodeHuff(AACSeq3(i).chl.sfc, 12, huffLUT);
            sfcr = decodeHuff(AACSeq3(i).chr.sfc, 12, huffLUT);

            Sl = decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook, huffLUT);
            Sr = decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook, huffLUT);

            frameFL = iAACquantizer(Sl', sfcl', AACSeq3(i).chl.G, frameType);
            frameFR = iAACquantizer(Sr', sfcr', AACSeq3(i).chr.G, frameType);

            frameFchl = iTNS(frameFL,frameType,AACSeq3(i).chl.TNScoeffs);
            frameFchr = iTNS(frameFR,frameType,AACSeq3(i).chr.TNScoeffs);

            frameF = [frameFchl frameFchr];
            frameT = iFilterbank(frameF, frameType, winType);

            x(s_start:s_end,1) = x(s_start:s_end,1) + frameT(:,1);
            x(s_start:s_end,2) = x(s_start:s_end,2) + frameT(:,2);
        end
    end

    %Remove any cliping (A small amount of samples is a little above 1 and below -1)
    x(x>1) = 1;
    x(x<-1) = -1;

    %Write the decoded waveform on the output file (48kHz sampling rate)
    audiowrite(fNameOut,x,48000);

end