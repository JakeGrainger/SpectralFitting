function [parameterArray,fvalArray,varParameterArray, invEH] = estimateParameter(obj,timeseries,Delta,fitIndex)
% estimateParameter  Estimate parameters of a spectral model from data.
%   estimateParameter(obj,timeseries,Delta,fitIndex) will estimate the
%   parameters of a model defined by obj such that they best describe the
%   timeseries provided. If timeseries is a matrix, then each column is
%   treated as a seperate timeseries, and a parameter estimate is provided.
%   Delta is the sampling interval of the timeseries and fitIndex is an
%   optional logical vector specifying which frequencies should be used in
%   the fitting.
%
%   For more information, see the Spectral Fitting Toolbox documentation.
arguments % error handling
    obj spectralFitting.SpectralModel
    timeseries {mustBeNumeric,mustBeNonempty,mustBeNonNan}
    Delta(1,1) {mustBeNumeric,mustBeNonempty,mustBeNonNan}
    fitIndex(:,1) = []
end
%% difference if applicable
if obj.useDifferencing
    timeseries = diff(timeseries,1,1); % differences the timeseries
end
% check fitting index is sensible
switch obj.objectiveMethod
    case 'bartlettLeastSquares'
        if isempty(fitIndex)
            fitIndex = (2:floor(obj.bartlettWindowSize/2))'; % if no fitIndex specified, then specify a fitIndex
        end
        testOmega = (1:floor(obj.bartlettWindowSize/2)+1)'; % set up a test array if fitIndex is specified to check it works
    otherwise
        if isempty(fitIndex)
            fitIndex = (2:floor(length(timeseries)/2))'; % if no fitIndex specified, then specify a fitIndex
        end
        testOmega = (1:floor(length(timeseries)/2)+1)'; % set up a test array if fitIndex is specified to check it works
end

try
    testOmega(fitIndex); % try fitIndex on testOmega
catch ME
    baseException = MException('estimateParameter:badfitsubscript',...
        'fitIndex is not compatible with provided timeseries and objectiveMethod.');
    baseException = addCause(baseException,ME);
    throw(baseException)
end
%% check compatibility of objective with fitting
fminconOpts = obj.fminconOpts;
switch obj.objectiveMethod
    case 'maximumLikelihood'
        if obj.useDifferencing % return error if differencing is tried with maximum likelihood
            error("%s not compatible with differencing.",obj.fitRoutine)
        end
end

%% set indexing
obj.d2IndexCon = reshape(obj.hessIndex,[],1); % used to specify which second derivatives are non-zero (probably should be a local variable and not a property)

%% start computations
switch obj.objectiveMethod
    case 'bartlettLeastSquares'
        omega = computeOmega(obj.bartlettWindowSize,Delta,fitIndex); % constuct a vector of angular frequencies, omega
        S = computeBartI(timeseries,Delta,fitIndex,obj.bartlettWindowSize); % get S, the spectral estimate
        objectiveFunction = @(parameter,record) LS(obj,S(:,record),omega,Delta,parameter); % make objective
    case 'maximumLikelihood'
        objectiveFunction = @(parameter,record) maximumLikelihood(obj,timeseries(:,record),Delta,parameter);
    case  'taperedDW'
        N = size(timeseries,1);
        filts=dpss(N,obj.taperBandWidth)'; hh1=filts(1,:)'; % dpss taper of bandwith 4
        HHK=zeros(N,1); % precomputation of h_t h_{t+\tau},described on page 16
        for mm=0:N-1
            HHK(mm+1)=sum(hh1(1:N-mm).*hh1(1+mm:N));
        end
        S = computeI(timeseries.*hh1*sqrt(size(timeseries,1)),Delta,fitIndex);
        objectiveFunction = @(parameter,record) taperedDW(obj,S(:,record),size(timeseries,1),Delta,fitIndex,parameter,HHK);
    otherwise
        % construct I
        omega = computeOmega(size(timeseries,1),Delta,fitIndex);
        S = computeI(timeseries,Delta,fitIndex);
        switch obj.objectiveMethod
            case 'leastSquares'
                objectiveFunction = @(parameter,record) LS(obj,S(:,record),omega,Delta,parameter);
            case 'Whittle'
                objectiveFunction = @(parameter,record) Whittle(obj,S(:,record),omega,Delta,parameter);
            case 'aliasedWhittle'
                objectiveFunction = @(parameter,record) aliasedWhittle(obj,S(:,record),omega,Delta,parameter);
            case 'debiasedWhittle'
                objectiveFunction = @(parameter,record) debiasedWhittle(obj,S(:,record),size(timeseries,1),Delta,fitIndex,parameter);
            otherwise
                ME = MException('estimateParameter:unknownObjectiveMethod',...
                    'Unknown objective method.');
                throw(ME)
        end
end
%% perform fitting
nRecord = size(timeseries,2);
parameterArray = nan(length(obj.parameter),nRecord);
fvalArray = nan(1,nRecord);
switch obj.fitRoutine
    case 'fminsearch'
        for iRecord = 1:nRecord
            boundedObjFun = @(parameter) objectiveBound(obj.lowerBound,...
                obj.upperBound, @(parameter) objectiveFunction(parameter,iRecord),parameter);
            [parameterArray(:,iRecord),fvalArray(:,iRecord)] =...
                fminsearch(boundedObjFun,obj.parameter,obj.fminsearchOpts);
        end
    case 'fmincon'
        lb = obj.lowerBound;
        ub = obj.upperBound;
        for iRecord = 1:nRecord
            [parameterArray(:,iRecord),fvalArray(:,iRecord)] =...
                fmincon(@(parameter) objectiveFunction(parameter,iRecord),...
                obj.parameter,[],[],[],[],lb,ub,[],fminconOpts);
        end
    case 'global'
        lb = obj.lowerBound;
        ub = obj.upperBound;
        for iRecord = 1:nRecord
            problem = createOptimProblem('fmincon',...
                'objective', @(parameter) objectiveFunction(parameter,iRecord),...
                'x0',obj.parameter,'lb',lb,'ub',ub,...
                'options',fminconOpts);
            gs = GlobalSearch('Display','off');
            [parameterArray(:,iRecord),fvalArray(:,iRecord)] = run(gs,problem);
        end
    case 'sequential'
        lb = obj.lowerBound;
        ub = obj.upperBound;
        problem = createOptimProblem('fmincon',...
            'objective', @(parameter) objectiveFunction(parameter,1),...
            'x0',obj.parameter,'lb',lb,'ub',ub,...
            'options',fminconOpts);
        gs = GlobalSearch('Display','off');
        [parameterArray(:,1),fvalArray(:,1)] = run(gs,problem);
        if nRecord > 1
            for iRecord = 2:nRecord
                [parameterArray(:,iRecord),fvalArray(:,iRecord)] =...
                    fmincon(@(parameter) objectiveFunction(parameter,iRecord),...
                    parameterArray(:,iRecord-1),[],[],[],[],lb,ub,[],fminconOpts);
            end
        end
    case 'debug'
        parameterArray = @(parameter) objectiveFunction(parameter,1);
    case 'oracleHess'
        parameterArray = obj.parameter;
end
% get uncertainty
if nargout > 2
    switch obj.objectiveMethod
        case 'debiasedWhittle'
            varParameterArray = nan(length(obj.parameter),length(obj.parameter),nRecord);
            invEH = nan(length(obj.parameter),length(obj.parameter),nRecord);
            for iRecord = 1:nRecord
                expHess = WhittleExpHess(obj,size(timeseries,1),Delta,fitIndex,...
                    parameterArray(:,iRecord));
                covLikelihood = covDebiasedWhittle(obj,size(timeseries,1),Delta,fitIndex,...
                    parameterArray(:,iRecord));
                varParameterArray(:,:,iRecord) = (expHess\covLikelihood)/expHess;
                invEH(:,:,iRecord) = inv(expHess);
            end
        otherwise
            varParameterArray = [];
            warning(['Variance not implemented for this objective function, '...
                'returning empty array.'])
    end
end

end

%% helper functions
function varargout = objectiveBound(lb,ub,objectiveFunction,parameter)
% function to put bounds on the parameter space, to be used by fminsearch
if any(parameter < lb) || any(parameter > ub)
    varargout = mat2cell(Inf(nargout,1),ones(nargout,1));
else
    [varargout{1:nargout}] = objectiveFunction(parameter);
end
end

function sqhess = fullHess(hess,p)
% function to convert to full square Hessian matrix
sqhess = zeros(p);
sqhess(tril(true(p))) = hess;
sqhess = sqhess + sqhess' - diag(diag(sqhess));
end

function fmix = crossMultiplyDerivative(deriv)
% function to cross multiply first derivative for use with upper triangular
% hessian
p = size(deriv,2);
f1 = repelem(deriv,1,p:-1:1);
f2 = repmat(deriv,1,p);
f2 = f2(:,reshape(tril(true(p)),p^2,1));
fmix = f1.*f2;
end
%% spectral estimates

function S = computeI(timeseries,Delta,fitIndex)
% function to calculate bartlett periodogram and corresponding omega
% S
S = Delta/(size(timeseries,1)*2*pi)*abs(fft(timeseries,[],1)).^2;
S = S(fitIndex,:);
end

function S = computeBartI(timeseries,Delta,fitIndex,windowSize)
% function to calculate bartlett periodogram and corresponding omega
% S
nSegments = floor(size(timeseries,1)/windowSize);
segmentPeriodograms = nan(windowSize,size(timeseries,2),nSegments);
for iSeg = 1:nSegments
    segmentPeriodograms(:,:,iSeg) = Delta/(windowSize*2*pi)*...
        abs(fft(timeseries((iSeg-1)*windowSize+1:iSeg*windowSize,:),[],1)).^2;
end
S = mean(segmentPeriodograms,3);
S = S(fitIndex,:);

end

function omega = computeOmega(N,Delta,index)
% function to compute omega
omega = 2*pi/Delta*(0:1/N:1/2);
omega = omega(index)';
end
%% objective functions

function [fval,deriv,hess] = LS(obj,S,omega,Delta,parameter)
% subfunction to calculate least squares objective
[density{1:nargout}] = parSpectralDensity(obj,omega,Delta,parameter); % get density and derivatives
diffS = density{1}-S; % calculate difference for use later
fval = sum(diffS.^2); % calculate least squares value
if nargout > 1
    deriv = 2*sum(density{2}.*diffS,1); % calculate derivative if required
end
if nargout > 2
    fmix = crossMultiplyDerivative(density{2});
    d2part = zeros(1,size(fmix,2));
    d2part(:,reshape(obj.hessIndex,[],1)) = sum(density{3}.*diffS,1);
    hess = fullHess(2*(d2part(:,reshape(tril(true(size(density{2},2))),[],1))+sum(fmix,1)),...
        size(density{2},2)); % calculate hessian if required
end
end

function [fval,deriv,hess] = generalWhittle(obj,density,S)
% subfunction to calculate whittle likelihood that takes a density argument
% so can be used by specific versions of the whittle likelihood
fval = sum(log(density{1}))+sum(S./density{1});
if nargout > 1
    deriv = sum(density{2}.*(1./density{1}-S./density{1}.^2),1);
end
if nargout > 2
    fmix = crossMultiplyDerivative(density{2});
    d2part = zeros(1,size(fmix,2));
    d2part(:,obj.d2IndexCon) = sum(density{3}.*(1./density{1}-S./density{1}.^2),1);
    hess = fullHess(d2part(:,obj.d2IndexUse)+...
        sum(fmix.*(-1./density{1}.^2+2*S./density{1}.^3),1),size(density{2},2));
end
end

function varargout = Whittle(obj,S,omega,Delta,parameter)
% subfunction to compute standard Whittle likelihood
[density{1:nargout}] = parSpectralDensity(obj,omega,Delta,parameter);
[varargout{1:nargout}] = generalWhittle(obj,density,S);
end

function varargout = aliasedWhittle(obj,S,omega,Delta,parameter)
% subfunction to compute aliased Whittle likelihood
[density{1:nargout}] = parAliasedSpectralDensity(obj,omega,Delta,parameter);
[varargout{1:nargout}] = generalWhittle(obj,density,S);
end

function varargout = debiasedWhittle(obj,S,N,Delta,index,parameter)
% subfunction to compute debiased Whittle likelihood
[density{1:nargout}] = parExpectedPeriodogram(obj,N,Delta,parameter,index);
[varargout{1:nargout}] = generalWhittle(obj,density,S);
end

function varargout = taperedDW(obj,S,N,Delta,index,parameter,taper)
% subfunction to compute tapered debiased Whittle likelihood
[density{1:nargout}] = parTaperedExpectedPeriodogram(obj,N,Delta,parameter,index,taper);
[varargout{1:nargout}] = generalWhittle(obj,density,S);
end

function varargout = maximumLikelihood(obj,timeseries,Delta,parameter)
% subfunction to calculate maximum likelihood for Gaussian process
[autoCov{1:nargout}] = parAutoCovariance(obj,length(timeseries),Delta,parameter);
cholSigma = chol(toeplitz(autoCov{1}));
invSigma = inv(cholSigma)/cholSigma';
varargout{1} = 0.5*(2*sum(log(diag(cholSigma))) + (timeseries)'*invSigma*(timeseries));
if nargout > 1
    E = invSigma * timeseries;
    F = (timeseries)' * invSigma;
    deriv = nan(1,size(autoCov{2},2));
    dSigma = cell(size(autoCov{2},2), 1);
    for iPar = 1:size(autoCov{2},2)
        dSigma{iPar} = toeplitz(autoCov{2}(:,iPar));
    end
    for iPar = 1:size(autoCov{2},2)
        deriv(1,iPar) = 0.5*(trace(invSigma * dSigma{iPar}) - F*dSigma{iPar}*E);
    end
    varargout{2} = deriv;
end
if nargout > 2
    hess = nan(size(autoCov{2},2));
    ind_mat = zeros(size(hess));
    ind_mat(tril(true(size(ind_mat)))) = 1:size(autoCov{3},2);
    for iPar = 1:size(autoCov{2},2)
        for jPar = 1:iPar
            h_part1 = dSigma{jPar}*invSigma*dSigma{iPar};
            index = ind_mat(iPar, jPar);
            d2Sigma = toeplitz(autoCov{3}(:,index));
            h_temp = 0.5*(-trace(invSigma*(h_part1-d2Sigma)) + F*(h_part1 - d2Sigma + dSigma{iPar}*invSigma*dSigma{jPar})*E);
            hess(iPar, jPar) = h_temp;
            hess(jPar, iPar) = h_temp;
        end
    end
    varargout{3} = hess;
end

end

%% expected hessians

function expHess = WhittleExpHess(obj,N,Delta,index,parameter)
% subfunction to calculate the expected hessian of the debiased Whittle
% likelhood
[density,derivative] = parExpectedPeriodogram(obj,N,Delta,parameter,index);
fmix = crossMultiplyDerivative(derivative);
expHess = fullHess(sum(fmix./(density).^2,1),size(derivative,2));
end

%% covariance of the debiased Whittle likelihood
function covLikelihood = covDebiasedWhittle(obj,N,Delta,index,parameter)
% function to calculate the covariance of the debiased Whittle likelihood
[EI,dEI] = parExpectedPeriodogram(obj,N,Delta,parameter,index);
a = -dEI./(EI.^2); % pxN
p = length(parameter);
covI = parCovPeriodogram(obj,N,Delta,parameter);
covIused = covI(index,index);
covLikelihood = nan(p,p);
for r = 1:p
    for s = r:p
        covLikelihood(r,s) = a(:,r)'*covIused*a(:,s);
        if r~=s
            covLikelihood(s,r) = covLikelihood(r,s);
        end
    end
end
end