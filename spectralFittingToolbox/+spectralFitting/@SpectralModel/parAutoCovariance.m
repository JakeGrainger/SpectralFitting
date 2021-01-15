function varargout = parAutoCovariance(obj,N,Delta,parameter)
% method to compute the autocovariance at lags 0,Delta,...,N-1*Delta for
% any parameter value
if isempty(obj.autoCovariance)
    M = 2*obj.chooseM(N); % times 2 to avoid periodicity in dft.
    omega = 2*pi/M/Delta*(0:floor(M/2))'; % corresponding Fourier frequencies
    [density{1:nargout}] = parAliasedSpectralDensity(obj,omega,Delta,parameter);
    varargout = cellfun(@(S) acf(S,N,Delta),density,"UniformOutput",false);
else
    [varargout{1:nargout}] = obj.autoCovariance(parameter,N,Delta);
end
end

function autocovariance = acf(S,N,Delta)
% subfunction that computes acf of columns of a matrix
S = [S;S(end-1:-1:2,:)];
oversampledAutocovariance = ifft(S,[],1)*2*pi/Delta;
autocovariance = real(oversampledAutocovariance(1:N,:));

end