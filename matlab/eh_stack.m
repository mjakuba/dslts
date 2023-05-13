classdef eh_stack < DslTsStack
    
    methods
        
        % constructor.
        function obj = eh_stack();
            obj@DslTsStack(1000,10);  % buflen, acclen
        end
        
        function m = apply_masks(obj,d,dbuf)            

            % hampel identifier (signed version).  This operates on
            % the buffer only.
                m = 0;
                m = m + obj.mask_hampel(d,dbuf,1e-5);

                % mask out positive values as uninteresting.
                % why does this break things?  Now it misses the
                % big peak?  Probably because the peak has lots of
                % positive values associated with it and so the
                % hampel identifier doesn't work properly.  This
                % should be in a separate stack anyway, or, since
                % this end up being used a lot, add the option
                % to ignore outliers of the wrong sign in the
                % hampel identifier.  No... that's unnatural
                % because the signal is not necessarily zero-mean.
                %m = m + obj.mask_abs_maximum(d,dbuf,0.0);
                
                
        end
        function f = apply_filters(obj,d,dbuf)
        % no filters, just fill with the output of the masking
        % stage.  Alternately we could simplify some accumulation
        % stages by doing a sort of pseudo-mask that actually
        % manipulates values.  So here we would make everything
        % that wasn't masked only as an outlier zero.  Kludgey.
        %f = obj.buffer.raw(obj.buffer.lst,1);
            
        %f = obj.filter_mask_not(d,dbuf,2^4);  % hampel bitmask is
        %2^4.
        % This is now done in the final stack.
            f = d(1);  
            
        end
        function accumulate(obj)
            
        % This is the unsigned version for now.  Really we just
        % want negative excursions.
        % obj.accumulate_tail_dt(100);
        % alter to do this in the final stack so we can ignore
        % non-negative there.
        end

    end
end
