function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
    
    xtrue = audioread(fNameIn); Ntrue = size(xtrue,1);
    x = padinput(xtrue,0);
    total_time = size(x,1)/48000;

    fprintf("Multimedia Systems:Part 3\n");

    AACSeq2 = AACoder2(fNameIn);
    
    fprintf("Compressed Encoding(V3):\n");
    tic; AACSeq3 = AACoder3(fNameIn, frameAACoded); toc;

    fprintf("Compressed Decoding(V3):\n");
    tic; xpred = iAACoder3(AACSeq3, fNameOut); toc;

    xtrue = x(1024:1024+Ntrue,:);
    xpred = xpred(1024:1024+Ntrue,:);

    err = xtrue - xpred;
    sig = xpred;

    SNR(1) = snr(sig(:,1), err(:,1));
    SNR(2) = snr(sig(:,2), err(:,2));

    sizeUncompressed = ByteSize(AACSeq2);
    sizeCompressed = ByteSize(AACSeq3);

    bitrateUncompressed = sizeUncompressed/total_time;
    sizeCompression = sizeUncompressed/sizeCompressed;

    bitrate = sizeCompressed/total_time;
    compression = bitrateUncompressed/bitrate;

    fprintf("Size Compression: size(AACSeq2)/size(AACSeq3) = %.2f\n",sizeCompression);
    fprintf("Bitrate Compression: bitrate(AACSeq2)/bitrate(AACSeq3) = %.2f\n",compression);
end


%Original Source: https://stackoverflow.com/questions/4845561/how-to-know-the-size-of-a-variable-in-matlab
function sizeBytes = ByteSize(in, fid)
    % BYTESIZE writes the memory usage of the provide variable to the given file
    % identifier. Output is written to screen if fid is 1, empty or not provided.
    
    if nargin == 1 || isempty(fid)
        fid = 1;
    end
    
    s = whos('in');
    %fprintf(fid,[Bytes2str(s.bytes) '\n']); 
    sizeBytes = s.bytes;
    end

    
    
    function str = Bytes2str(NumBytes)
    % BYTES2STR Private function to take integer bytes and convert it to
    % scale-appropriate size.
    
    scale = floor(log(NumBytes)/log(1024));
    switch scale
        case 0
            str = [sprintf('%.0f',NumBytes) ' b'];
        case 1
            str = [sprintf('%.2f',NumBytes/(1024)) ' kb'];
        case 2
            str = [sprintf('%.2f',NumBytes/(1024^2)) ' Mb'];
        case 3
            str = [sprintf('%.2f',NumBytes/(1024^3)) ' Gb'];
        case 4
            str = [sprintf('%.2f',NumBytes/(1024^4)) ' Tb'];
        case -inf
            % Size occasionally returned as zero (eg some Java objects).
            str = 'Not Available';
        otherwise
           str = 'Over a petabyte!!!';
    end    
end