function varargout = parTaperedExpectedPeriodogram(obj,N,Delta,parameter,index,taper)
% methods to compute the expected periodogram, for any parameter value
[autoCov{1:nargout}] = parAutoCovariance(obj,N,Delta,parameter);
varargout = cellfun(@(C) expTI(C,Delta,index,taper),autoCov,"UniformOutput",false);

end

function EI = expTI(autocovariance,Delta,index,taper)
% subfunction to calculate expected periodogram for each column of a
% matrix
EI = Delta/(2*pi)*...
    real(2*fft(autocovariance.*taper,[],1)-autocovariance(1,:));
EI = EI(index,:); % include desired frequencies
end