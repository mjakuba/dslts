classdef ehanom_stack < DslTsStack
    
    methods
        
        % constructor.
        function obj = ehanom_stack();
            obj@DslTsStack(2,10);  % buflen, acclen
        end
        
        function m = apply_masks(obj,d,dbuf)            

        % Ignore positive anomalies.
            m = obj.mask_abs_maximum(d,dbuf,0.0);
                
        end
        function f = apply_filters(obj,d,dbuf)
        % no filters, just fill with the output of the masking
        % stage.  Alternately we could simplify some accumulation
        % stages by doing a sort of pseudo-mask that actually
        % manipulates values.  So here we would make everything
        % that wasn't masked only as an outlier zero.  Kludgey.
        %f = obj.buffer.raw(obj.buffer.lst,1);
            
            f = obj.filter_mask_not(d,dbuf,2^4);  % hampel bitmask is 2^4.
            
        end
        function accumulate(obj)
          
        % This is the unsigned version for now.  Really we just
        % want negative excursions.  Fixed this by using an
        % additional stack and mask.
            obj.accumulate_tail_dt(100);
        end

    end
end
