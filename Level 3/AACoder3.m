
function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
    %%Read the audio from the input file and padd it accordingly
    x = padinput(audioread(fNameIn),1);
    N = max(size(x));
    huffLUT = loadLUT();
    
    AACSeq2 = AACoder2(fNameIn);
    K = length(AACSeq2);
    
    for i=1:K-1
        AACSeq3(i).frameType = AACSeq2(i).frameType;
        AACSeq3(i).winType = AACSeq2(i).winType;
    
        AACSeq3(i).chl.TNScoeffs = AACSeq2(i).chl.TNScoeffs;
        AACSeq3(i).chr.TNScoeffs = AACSeq2(i).chr.TNScoeffs;
    end
    
    frames = []; K = 0;
    inds = 1-1024:1024;
    
    while not( inds(end) == N )
        inds = inds + 1024;
        K = K+1;

        frames = [frames x(inds,:);];
    end
    
    
    for i=1:K-1
        frame = frames(:,[i*2-1 i*2]);
        switch i
        case 2
            frameTprev2 = zeros(2048,2);
            frameTprev1 = frames(:,[(i-1)*2-1 (i-1)*2]);
        case 1
            frameTprev2 = zeros(2048,2);
            frameTprev1 = zeros(2048,2);
        otherwise
            frameTprev2 = frames(:,[(i-2)*2-1 (i-2)*2]);
            frameTprev1 = frames(:,[(i-1)*2-1 (i-1)*2]);
        end
    
        frameType = AACSeq3(i).frameType;
        SMRleft = psycho(frame(:,1),frameType,frameTprev1(:,1),frameTprev2(:,1));
        SMRright = psycho(frame(:,2),frameType,frameTprev1(:,2),frameTprev2(:,2));
        
        [Sl, sfcl, AACSeq3(i).chl.G] = AACquantizer(AACSeq2(i).chl.frameF, frameType, SMRleft);
        [Sr, sfcr, AACSeq3(i).chr.G] = AACquantizer(AACSeq2(i).chr.frameF, frameType, SMRright);
    
        if strcmp(frameType,"ESH")
            [AACSeq3(i).chl.sfc, ~] = encodeHuff(sfcl(:), huffLUT, 12);
            [AACSeq3(i).chr.sfc, ~] = encodeHuff(sfcr(:), huffLUT, 12);
        else
            [AACSeq3(i).chl.sfc, ~] = encodeHuff(sfcl(:), huffLUT, 12);
            [AACSeq3(i).chr.sfc, ~] = encodeHuff(sfcr(:), huffLUT, 12);
        end
    
        [AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook] = encodeHuff(Sl(:), huffLUT );
        [AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook] = encodeHuff(Sr(:), huffLUT );
    
    end

    save(fnameAACoded, 'AACSeq3');
end
    

