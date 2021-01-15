function varargout = parExpectedPeriodogram(obj,N,Delta,parameter,index)
% methods to compute the expected periodogram, for any parameter value
[autoCov{1:nargout}] = parAutoCovariance(obj,N,Delta,parameter);
varargout = cellfun(@(C) expI(C,N,Delta,index),autoCov,"UniformOutput",false);

end

function EI = expI(autocovariance,N,Delta,index)
% subfunction to calculate expected periodogram for each column of a
% matrix
EI = Delta/(2*pi)*...
    real(2*fft(autocovariance.*(1-(0:N-1)/N)',[],1)-autocovariance(1,:));
EI = EI(index,:); % include desired frequencies
end