function covI = parCovPeriodogram(obj,N,Delta,parameter)
% method to calculate the covariance of the periodogram
M = obj.chooseM((2*N-2)+1); % because we need (2N-2)+1 points, and may need over sampling
omega = pi/Delta/M*[0:1:(M-1), -M:1:-1]'; % 2M, avoids possible issues with not having enough lags
density = obj.parAliasedSpectralDensity(omega,Delta,parameter);
q = density.*exp(1i*Delta*(N-1)*omega);
Q_long = fft(q)*pi/M/Delta; % is P 2 pi over 2 M
Q_vec = Q_long((0:2*N-2)+1); % recover desired lags
Q = hankel(Q_vec(1:N),Q_vec(N:end)); % makes hankel matrix
covI = abs((Delta/2/pi/N)*N^2*ifft2(Q)).^2;
end