classdef dehdt_stack < DslTsStack
    
    methods
        
        % constructor.
        function obj = dehdt_stack();
            obj@DslTsStack(20,0);  % buflen, acclen
        end
        
        function m = apply_masks(obj,d,buf)
            
            m = 0;
            b = obj.mask_depth_min(d,buf,1000);  % For compiled classes this
                                           % can be part of the constructor.
            m = m+b;
            
        end
        function f = apply_filters(obj,d,dbuf)

        % for Eh processing we need dEhdt - get this from the
        % filter.
        %f = obj.filter_diff(d,dbuf);
        % straight diff causes implosion of the scale estimate in
        % the later hampel identifier because of quantization in
        % the raw signal.
            f = obj.filter_diff_poly(d,dbuf,1,10);
        end
        
        function accumulate(obj)
        
        % do nothing.  Or alternately we could do a tail here,
        % without the intermediate step of a hampel identifier.
        % This would capture the best values so far, without 
        % attempting at all to define a threshold.  We probably
        % want to see whether we can do a batch hampel identifier
        % as well (by estimating the median and maddev and
        % generating a threshold that is then applied to incoming
        % data (and accumulated data in the tail that fails the
        % (evolving) threshold is discarded.
            
        end
            
    end
end
