function SNR = demoAAC2(fNameIn, fNameOut)

    xtrue = audioread(fNameIn); 
    trueSize = size(xtrue,1);

%%Encoding
    disp("Encoding starting"); tic();
    AACSeq2 = AACoder2(fNameIn); toc();

%%Decoding  (In xpred we remove the padded pre and post zeroes)
    disp("Decoding starting"); tic();
    xpred = iAACoder2(AACSeq2, fNameOut); toc();    xpred = xpred(1025:end,:);  xpred = xpred(1:trueSize,:);

%%Calculate SNR
    err = xtrue - xpred;
    sig = xtrue;

    SNR(1) = snr(sig(:,1),err(:,2));
    SNR(2) = snr(sig(:,2), err(:,2));

end
