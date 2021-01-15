function varargout = parAliasedSpectralDensity(obj,omega,Delta,parameter)
% function to compute the aliased spectral density, will approximate if not
% specified in object property, for any parameter value.
if isempty(obj.aliasedSpectralDensity) % aliasedSpectralDensity not provided
    K = obj.chooseK(Delta,parameter);
    % construction of aliasedSDF
    kArray = 2*pi/Delta*(-K:K); % make array of k values to add on to omega
    [varargout{1:nargout}] = obj.parSpectralDensity(omega+kArray(1),...
        Delta,parameter); % compute at first k value
    for k = 2:length(kArray)
        [kDensity{1:nargout}] = obj.parSpectralDensity(omega+kArray(k),...
            Delta,parameter); % temporarily store density at k
        for iOut = 1:nargout
            varargout{iOut} = varargout{iOut}+...
                kDensity{iOut}; % add to running total
        end
    end % loop may seem inefficient, but the memory cost of vectorising is worse
else % aliasedSpectralDensity has been provided
    [varargout{1:nargout}] = obj.aliasedSpectralDensity(parameter,omega,Delta);
end

end