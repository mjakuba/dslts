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

classdef obsprofile_stack < DslTsStack
    
    methods
        
        % constructor.
        function obj = obsprofile_stack();
            obj@DslTsStack(1,200);  % buflen, acclen.  Accumulator
                                   % is a trend estimator and needs
                                   % only two outputs.  Does it
                                   % need intermediate variables?
                                   % Almost certainly.
        end
        
        function m = apply_masks(obj,d,dbuf)            

        % cheating here by using a known-good depth band for
        % estimating background profile above NBP on abe126.
        % jakuba07 implements a somewhat successful automatic means
        % of segmenting some data above the NBP from the descent.
            m = 0;
            m = m + obj.mask_depth_minmax(d,dbuf,2120,2220);
            
            % also ignore data after the descent is completed.
            m = m + obj.mask_time_max(d,dbuf,12500);  % relative time.
            
            
        end
        function f = apply_filters(obj,d,dbuf)

        % no filtering, just ignore masked values.
            f = obj.filter_mask_not(d,dbuf,0);
            
        end
        function accumulate(obj)

        % iterative computation of trendline from unmasked data.
        % working here.
        % either do recursive LS fitting of line or something
        % completel hackish.  Not sure how to get noise estimate
        % needed for prediction interval using RLS but should be
        % straightforward if posed as a recursive estimation
        % problem.
        % key discovery for interface today is accumulators will be
        % their own class or similar.  In C++ we will use
        % boost:accumulator for this.
            obj.accumulate_lsline();
                        
        end

    end
end
