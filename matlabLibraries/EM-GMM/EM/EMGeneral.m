function [expectations theta logLike] = EMGeneral(X,expectations,theta,handles,options)

options = setDefaultOptions(options);

if ~isstruct(theta) && ~isstruct(expectations)
	error('Neither latent variables nor parameters initialized.')
end

% if expectations is initialized, calculate intial parameters 
if isstruct(expectations)

end	

logLike = 0;

for iter = 1:options.numIters
	% E STEP
	expectations = handles.estep(X,theta);

	% M STEP
	[theta logLike(iter+1)] = handles.mstep(X,expectations,theta);
	if options.verbose > 0
		fprintf('%d: %.4f\n',iter,logLike(iter+1));
	end

    if abs((logLike(iter+1)-logLike(iter))/logLike(iter+1))<10^-6
		fprintf('Optimization converged at iteration %d.\n',iter);
        break
    end

	if options.plotting
		handles.plotting(X,expectations,theta);
	end
end

end % end function


function options = setDefaultOptions(options)

if ~isfield(options,'verbose')
	options.verbose = 1;
end

if ~isfield(options,'numIters')
	options.numIters = 1000;
end

if ~isfield(options,'saveHist')
	options.saveHist = 0;
end

if ~isfield(options,'plotting')
	options.plotting = 0;
end

end % end function
