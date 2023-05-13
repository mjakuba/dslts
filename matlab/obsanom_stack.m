% anomaly relative to background profile.  profile is computed
% elsewhere (in a previous stack).

classdef obsanom_stack < DslTsStack
    
    methods
        
        % constructor.
        function obj = obsanom_stack();
            obj@DslTsStack(10,200); 
        end
        
        function m = apply_masks(obj,d,dbuf)            

        % no masks to apply.
            m = 0;
            
        end
        function f = apply_filters(obj,d,dbuf)

        % detrend.
            f = obj.filter_detrend_cov(d,dbuf);
            
        end
        function accumulate(obj)

        % could accumulate representative sample of NBP data, or
        % most intense NBP data, or do hampel identifier for BP
        % data.
        % for now, nothing.
        % or could collect (widest) intervals of nbp detection
        % that sort of data is of limited interest in the near
        % term, so implement the BP detector here.  In jakuba07
        % this requires that the data be both outside the
        % prediction interval (to be assessed as NBP) and
        % identified as an outlier by a hampel identifier.
                        
        end

    end
end
