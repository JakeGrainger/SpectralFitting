classdef SpectralModel
    % Model containing parametric forms for the spectral density function, alongside parameters and bounds for the parameter space
    % 
    properties
        % stores handle for the spectral density function
        spectralDensity function_handle
        % stores handle for the auto covariance
        autoCovariance function_handle
        % stores handle for the aliased sdf
        aliasedSpectralDensity function_handle
        % lower bounds for parameter space
        lowerBound(:,1) {mustBeNumeric}
        % upper bounds for parameter space
        upperBound(:,1) {mustBeNumeric}
        % parameter vector
        parameter(:,1) {mustBeNumeric}
        % index of elements of the upper triangular hessian that may be
        % non-zero
        hessIndex {mustBeNumericOrLogical}
    end
    properties (Access = public)
        % minimum number of bins used in integrals
        nIntegralBins(1,1) {mustBeNumeric} = 4096
        % aliasing tolerance
        aliasTolerance(1,1) {mustBeNumeric} = 1e-6
        % objective function for fitting
        objectiveMethod {mustBeMember(objectiveMethod,...
            {'debiasedWhittle', 'Whittle', 'aliasedWhittle',...
            'leastSquares', 'bartlettLeastSquares', 'maximumLikelihood', 'taperedDW'})} = 'debiasedWhittle'
        % indicator for if differencing should be used
        useDifferencing(1,1) logical {mustBeNumericOrLogical} = false
        % fitting routine to be used
        fitRoutine {mustBeMember(fitRoutine,{'fminsearch','fmincon',...
            'global','sequential', 'debug', 'oracleHess'})} = 'fminsearch'
        % bartlett smoothing window
        bartlettWindowSize(1,1) {mustBeNumeric} = 128
        % taper band width
        taperBandWidth = 4
        
    end
    properties (Access = private)
        % index for constructing d2
        d2IndexCon
        % index for using d2
        d2IndexUse
    end
    properties
        % options for fminsearch
        fminsearchOpts = optimset('GradObj','on','MaxFunEvals',100000,...
            'MaxIter',100000,'TolFun',1e-8,'TolX',1e-8)
        fminconOpts = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
            'SpecifyObjectiveGradient',true,'HessianFcn','objective',...
            'MaxIterations',100000,'Display','off')
    end
    methods
        function obj = SpectralModel(spectralDensity,parameter,...
                lowerBound,upperBound)
            % Constructor function that takes as input a spectralDensity,parameter and optional arguments lowerBound and upperBound
            arguments
                spectralDensity function_handle
                parameter(:,1) {mustBeNumeric}
                lowerBound(:,1) = -Inf(size(parameter))
                upperBound(:,1) = Inf(size(parameter))
            end
            if nargin < 3
                warning('Bounds not specified. Setting bounds to [-Inf,Inf] for each parameter.')
            elseif nargin < 4
                warning('Upper bound not specified. Setting to Inf for each parameter.')
            end
            obj.spectralDensity = spectralDensity;
            obj.parameter = parameter;
            obj.lowerBound = lowerBound;
            obj.upperBound = upperBound;
            obj.hessIndex = tril(true(length(parameter))); % defaults to calcualating the whole Hessian matrix (without assuming zeros)
            obj.d2IndexCon = reshape(obj.hessIndex,[],1); % gets set in estimation routine anyway
            obj.d2IndexUse = reshape(tril(true(length(parameter))),[],1); % used for indexing in optimisation, only depends on size of parameter vector
        end
        varargout = computeSpectralDensity(obj,omega,Delta)
        varargout = computeAliasedSpectralDensity(obj,omega,Delta)
        varargout = computeAutoCovariance(obj,N,Delta)
        varargout = computeExpectedPeriodogram(obj,N,Delta)
        varargout = simulateGaussianProcess(obj,N,Delta);
        varargout = estimateParameter(obj,timeseries,Delta);
        series    = simulateTucker(obj,N,Delta,M,useRandAmplitudes);
        covI      = computeCovPeriodogram(obj,N,Delta,parameter);
        corrI     = computeCorrPeriodogram(obj,N,Delta,parameter);
    end
    methods % plotting
        varargout = plotSdf(obj,N,Delta,varargin)
        varargout = plotdBSdf(obj,N,Delta,varargin)
        varargout = plotAliasedSdf(obj,N,Delta,varargin)
        varargout = plotdBAliasedSdf(obj,N,Delta,varargin)
        varargout = plotAcv(obj,N,Delta,varargin)
        varargout = plotEI(obj,N,Delta,varargin)
        varargout = plotdBEI(obj,N,Delta,varargin)
        varargout = plotCovI(obj,N,Delta,varargin)
        varargout = plotCorrI(obj,N,Delta,varargin)
    end
    methods (Access = private) % contruction of useful parameters
        function M = chooseM(obj,N)
            % compute M based on N and preferences in the object
            M = N*ceil(obj.nIntegralBins/N);
        end
        function K = chooseK(obj,Delta,parameter)
            % compute K based on Delta and preferences in the object
            K = 1;
            while parSpectralDensity(obj,(2*K+1)*pi/Delta,Delta,parameter) > obj.aliasTolerance
                K = K + 1;
            end
        end
    end
    methods (Access = public) % for optimisation
        varargout = parSpectralDensity(obj,omega,Delta,parameter)
        varargout = parAliasedSpectralDensity(obj,omega,Delta,parameter)
        varargout = parAutoCovariance(obj,N,Delta,parameter)
        varargout = parExpectedPeriodogram(obj,N,Delta,parameter,index)
        covI      = parCovPeriodogram(obj,N,Delta,parameter);
        varargout = parTaperedExpectedPeriodogram(obj,N,Delta,parameter,index,taper)
    end
end