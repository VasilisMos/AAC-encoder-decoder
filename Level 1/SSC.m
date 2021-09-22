function frameType = SSC(frameT, nextFrameT, prevFrameType)

    switch prevFrameType
    case "LSS"
        frameType = "ESH";
    case "LPS"
        frameType = "OLS";
    otherwise %ESH || OLS
        nextLeft = nextFrameT(:,1);
        nextRight = nextFrameT(:,2);

        nextLeftType = checkNextFrame(nextLeft);
        nextRightType = checkNextFrame(nextRight);

        nextType = decideType(nextLeftType,nextRightType);

        if strcmp(prevFrameType,"OLS")
            if strcmp(nextType,"OLS")
                frameType = "OLS"; return;
            else
                frameType = "LSS"; return;
            end
        end
        if strcmp(prevFrameType,"ESH")     
            if strcmp(nextType,"OLS")
                frameType = "LPS"; return;
            else
                frameType = "ESH"; return;
            end
        end
    end
end

%Implements the 
function [nextFrameType] = checkNextFrame(nextFrameT)

    %Filter coefficients
    n = [0.7548 -0.7548];
    dn = [1 -0.5095];

    indT = 576:703;
    l = [];

    for i=1:8
        l = [ l; filter(n,dn,nextFrameT(indT))' ];
        indT = indT + 128;
    end

    sl2 = sum(l.^2,2);
    dsl = zeros(length(sl2),1);

    for i =2:length(sl2)
        dsl(i) = sl2(i)/( (1/(i-1)) * sum(sl2(1:i-1)) );
    end

    for i =2:length(sl2)
        if sl2(i) > 0.001 && dsl(i) > 10
           nextFrameType = "ESH";
           return; %There is no need to search other subframes
        end
    end

    nextFrameType = "OLS"; %Not a single subrame meets the "ESH" requirements
end

%%"Hard Code" of the lookup Table given on the assignment pds (Page 5)
function [frameType] = decideType(frameType1,frameType2)
    switch frameType1
    case "OLS"
        switch frameType2
        case "OLS"
            frameType = "OLS";
        case "LSS"
            frameType = "LSS";
        case "LPS"
            frameType = "LPS";
        case "ESH"
            frameType = "ESH";
        end
    case "LSS"
        switch frameType2
        case "OLS"
            frameType = "LSS";
        case "LSS"
            frameType = "LSS";
        case "LPS"
            frameType = "ESH";
        case "ESH"
            frameType = "ESH";
        end
    case "LPS"
        switch frameType2
        case "OLS"
            frameType = "LPS";
        case "LSS"
            frameType = "ESH";
        case "LPS"
            frameType = "LPS";
        case "ESH"
            frameType = "ESH";
        end
    case "ESH"
        frameType = "ESH";
    otherwise
        frameType = "ERROR";
    end
end
