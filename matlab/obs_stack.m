% potential density.
%
% Used as a covariate for other scalars.  
%
% Depends on S,T,P, and latitude.
%
% So there are multiple inputs.  Makes no sense to put non-general
% processing of data here.  So this should operate on downstream
% stuff with sigma as a covariate.  Often will need to interpolate
% the output of one chain into another when the output of one is a
% covariate for another.  Probably need a utility
% function for doing that.  But how?  Should the dependent chain be
% delayed somehow?  Or we could ZOH for a maximum interval (or
% indefinitely), just assuming that the stacks run close to
% synchronously.  That's the right strategy for now, given limited time.

classdef obs_stack < DslTsStack
    
    methods
        
        % constructor.
        function obj = obs_stack();
            obj@DslTsStack(20,2);  % buflen, acclen.  Accumulator
                                   % is a trend estimator and needs
                                   % only two outputs.  Does it
                                   % need intermediate variables?
                                   % Almost certainly.
        end
        
        function m = apply_masks(obj,d,dbuf)            

            m = obj.mask_depth_min(d,dbuf,100); % ignore surface data.
            
        end
        function f = apply_filters(obj,d,dbuf)

        % De-spike
            f = obj.filter_median(d,dbuf,10);  % 10-pt per
                                               % jakuba07.  Seems a
                                               % lot?
            
        end
        function accumulate(obj)

        % do nothing here.  trend line computation has to happen in
        % a separate stack because we want to mask that data by
        % depth interval.
        end

    end
end
