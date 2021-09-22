function frameFout = iTNS(frameFin, frameType, TNScoeffs)
    
    switch frameType
    case "ESH"
        for i=1:8
            frameFout(:,i) = filter([1],[1 -(TNScoeffs(:,i)')],frameFin(:,i));
        end
    otherwise
        frameFout = filter([1],[1 -(TNScoeffs')],frameFin);
    end
end