% TODO (relative to prototype here)
% * filtered value and masked value should be part of main buffer.
%   No need for separate buffers and associated housekeeping.
% * Passing of processing stack and associated parameters.  Maybe
%   just ignore this for now because the syntax is language
%   dependent anyway.

% boost accumulator is chained -- determines its own dependencies
% but this chaining is fixed and dependent on the particular
% output desired and not on input parameters.  Might be a model for chaining
% masks and filters, and perhaps whole stacks, though no plans to
% allow mixing the order of masks, filters, and accumulators with
% stacks.  Except that filters and masks don't depend on one
% another, we just often want to use several.

% some global variables
%global FILTER_DIFF; FILTER_DIFF=0;
%global MASK_ABS_MAX; MASK_ABS_MAX=0;

% class definition.  handle class to make property access easier.
classdef DslTsStack < handle 

        properties
            
            % initialization data.
            % could conceivably make a default class that took as
            % initialization parameters IDs to the functions to
            % apply to data.  This might avoid having to compile
            % classes for the most basic processing.  Would not
            % support dependent instances of processing chains but
            % those could always be constructed as compiled code.
            % Actually that might still be supported because
            % compiled code will always be required in programs
            % that use this class.
            buflen
            acclen 
            % @@@ acclen is not general.  Some accumulators have
            % the number of samples to keep as a sole parameter.
            % Others have many more parameters, and may not keep
            % any samples.  Accumulators are really better
            % expressed as a class than as functions.  So the
            % method accumulate would add data to an existing
            % accumulator.  There would be some wrapper method to
            % pull the associated result out.  The result could be
            % some samples, or might be a mean, or parameters of a
            % trend line, etc.
            
            % data in circular buffer
            buffer
            
            % accompanying buffer of mask (idenfied via bitmask;
            % zero indicates not masked).
            %mask
            % acommpanying buffer of filtered values (primary data
            % only, covariates are not filtered).
            %filtered
            % no - combine mask and filtered value into the main
            % buffer to avoid indexing issues.
            
            % Variable-size output of the accumulation step.
            % How are we going to initialize this?  Just make it a
            % fixed initialization parameter to the class constructor.
            accumulated 
            num_accumulated
            num_samples 
            trend_intercept
            trend_slope
            pred_interval_min
            pred_interval_max
            
        end

        methods (Abstract)
            
            % available masks are always the same, but not
            % all are applied.  So use a bitmask when instantiating
            % the class to set which are applied.  Then we need
            % the parameters for each.  Do these at construction
            % as well.  Or we can define a bunch of useful ones
            % elsewhere and let subclasses use them.
            apply_masks(obj,d,buf)
            
            % available filters are always the same, but the order
            % of application varies and not all are applied 
            % (generally few will be applied).
            % More awkward to specify these because order is also
            % 
            % Right - so let subclass define this method using 
            % the fully defined methods in this superclass.
            apply_filters(obj,d,buf)
            % This is "pre-processing" and only 
            
            % Incremental statistics... What do we call this?  
            % Adding to accumulators
            % long-term statistics
            accumulate(obj)
            
            % Combined statistics
            % short-term statistics that require long-term stats
            % to compute, e.g. detrended data.  Output is a
            % processed version of the buffer.
            % combine(obj)
            % wouldn't this just be a filter?  Except that we need
            % the long term statistics to compute it, so it is
            % dependent on the accumulation step.
            %filter2(obj)
            
            % Second accumulation?  E.g. tail of detrended data?
            %accumulate2(obj);
            
            % ... may want to run through this thing again.  If so,
            % then shouldn't that require two instantiations of
            % this class, the output of one feeding the other?
            %
            % e.g. Implementation for outlier ID on de-trended OBS
            % data:
            % 1) pre-processing up through line-fitting in the
            % accumulate step
            % 2) Feed filtered data from (1) along with filter parameters
            % for linear-detrending in filter stage of (2), then
            % windowed hampel identification and finally collect
            % the tail.  Difficulty here is that linear-detrending
            % parameters evolve with time.  This will corrupt the
            % long-term statistics, but that is what it is.  So
            % parameters have to go in as floating-point
            % covariates?
            %
            % Not clear how to get even static parameters in here
            % for given filters, masks, etc, unless new classes
            % have to be derived for every processing chain.  Want
            % to pass something like 
            % [maskY,p1],[filterX,p1,p2,p3],[accumulatorB,p1,p2]
            % to the constructor.  Or, since it will be common to
            % want to apply several masks and filters in order, to
            % pass a list.  Then why not specify the entire
            % processing chain this way and not limit the number of
            % masks and filters applied?  (Accumulator output
            % doesn't have the same dimensionality so it probably
            % cannot appear throughout).  Also limit this to all
            % masks being applied before all filters, and all
            % filters before accumulation to limit implementation
            % complexity.
            %
            % This is probably pretty easy to implement in matlab,
            % but probably difficult though likely possible in C++.
            %
            % the de-trender could just rely on a particular
            % labling for the covariates, and the class could be
            % informed of what to expect as covariates at
            % construction (and fail if it doesn't see what it
            % expects, based on what filters it has been told it
            % will apply).
            %
            % Use compile-time defined derived classes to start and
            % verify the concept.  Then work on interface.
            
        end
            
        % This maybe is where we define all the useful filters
        % etc. that derived classes devoted to individual 
        % tracers can use.
        %
        % Also define the fundamental input and output.
        % Input: data into buffer along with required and optional
        % covariates
        %
        % Output: output from various processing stages as vectors
        % and strings.
        methods
            
            % constructor
            % function obj = DslTsStack(buflen,masks_to_apply,...);
            function obj = DslTsStack(buflen,acclen);
                
                datalen = 8;  % maybe this is the place to define
                              % the number of floating point covariates.
                obj.buflen = buflen;
                obj.acclen = acclen;

                % scalar plus covariates
                % @@@ covariates?
                % required: time, depth, northing, easting
                % optional: potential density, survey stage or
                % (multiple) other integer identifiers.
                % * integer identifiers, perhaps just one, that can
                % be used by custom filters, or just make this a
                % binary input that can be used to directly mask
                % data, e.g. to only pay attention to a single
                % survey stage.  Multiple indicators could be
                % combined external to the chain.
                % * altitude -- inevitably going to want this.
                % * one additional floating-point covariate.  Do we
                % ever want more?
                obj.buffer = circVBuf(int64(buflen),int64(datalen+1));
                % +1 for filtered value.  Mask is a bitmask so gets
                % combined with input mask.
                
                % pre-mask buffer for unfilled values.
                d0 = [NaN NaN NaN NaN NaN NaN 2^32 NaN NaN]; 
                for n=1:buflen
                    obj.buffer.append(d0);
                end
                
                % accompanying buffer of mask (idenfied via bitmask;
                % zero indicates not masked).
                %obj.mask = circVBuf(int64(buflen),int64(1));
                % Initialize mask to indicate all empty values
                %obj.mask.append(ones(buflen,1));
                % normally would have to reset location of
                % pointer.  cheating here.
                % probably make the mask part of the main buffer.
                % No point even allowing the possibility of getting
                % them mixed up.
                %
                % acommpanying buffer of filtered values (primary data
                % only, covariates are not filtered).
                % Only 1 filtered output per instance of this
                % class.  If intermediate or other products are
                % desired, instantiate another class.
                %obj.filtered = circVBuf(int64(buflen),int64(1));
                
                % in final implementation accumulators will
                % necessarily be their own classes.  HAcking this
                % for now by adding all data and parameters needed
                % by the union of all implemented accumulators.
                obj.accumulated = NaN*ones(acclen,datalen+1);  % +1
                                                               % for filterd
                obj.num_accumulated = 0;
                obj.trend_intercept = NaN;
                obj.trend_slope = NaN;
                obj.pred_interval_min = NaN;
                obj.pred_interval_max = NaN;
                obj.num_samples = 0;  % total number of samples
                                      % processed by accumulator.
            end
            
            % Consume data
            function add_data(obj,data)
               
                assert(length(data)==8);
                % assert(sum(isnan(data)) == 0);  Will need NaNs
                % for e.g. unused covariate.
                % scalar,time,northing,easting,depth,altitude,external_mask,scalar_covariate)
                % ?more (or perhaps unlimited) scalar covariates?
                
                % extend data and initialize new values.  Probably
                % data should be a structure.  Yes definitely a
                % structure or its own class, but can't use that
                % with circVBuf so just use a vector for now.
                d(1:8) = data;
                d(9) = data(1); % Default is no filtering in which
                              % case this is just the input.
                
                %disp('append data');
                %obj.buffer.append(buf(:)');
                
                % process new data....
                
                % masks.  
                % Always apply the input mask.  This will always be
                % done first (it could be done with the output of
                % one instance feeding another if required).
                m0 = d(7);
                dbuf = obj.buffer.raw(obj.buffer.fst:obj.buffer.lst,:);
                ma = apply_masks(obj,d,dbuf);
                m = m0+ma;  % bitmask.
                d(7) = m;

                % add masks to buffer
                %obj.mask.append(m);
                
                % filters
                f = apply_filters(obj,d,dbuf);
                % add filtered value to buffer
                %obj.filtered.append(f);
                d(9) = f;
                
                % append data
                obj.buffer.append(d);
                
                % accumulate data (iterative statistics and
                % estimators)
                % WTF? Cannot seem to make accumulated persist
                % between calls to add_data.  Other properties do.
                accumulate(obj);

                
            end
            
            % output data from various stages. 
            
            % raw contents of buffer
            function print_buffer(obj,n)
                disp('raw contents of circular buffer:');
                obj.buffer.raw
            end
            
            % Statistical processing functions.
            
            % Example absolute mask - data below a minimum is
            % masked
            % reserve 2^0 for input mask.
            % reserve 2^32 for unfilled buffer.
            function b = mask_abs_minimum(obj,d,buf,MIN)
                b = 0;
                if d(1) < MIN
                    b = 2^1;
                end
            end
            function b = mask_abs_maximum(obj,d,buf,MAX)
                b = 0;
                if d(1) > MAX
                    b = 2^2;
                end
            end
            function b = mask_depth_min(obj,d,buf,MIN);
                b = 0;
                if d(5) < MIN
                    b = 2^3;
                end
            end
            function b = mask_hampel(obj,d,buf,PF)
            % is the newest data an outlier?  Probably should use 
            % data already identified as outlier, but ignore all
            % other masks.
                
                b=0;
                %ii = ~logical(buf(:,7));  % index into unmaksed values
                ii = buf(:,7) ~= 2^4;  % already marked outliers
                                       % should be included.
                
                % ignore data until there are enough values;
                if sum(ii) < obj.buflen/2 % arbitrary threshold
                    disp('ignoring values buffer insufficiently full');
                    return;
                end
                
                % ignore data if latest is masked.
                if d(7) ~= 0
                    disp('latest data is masked');
                    return
                end
                
                % need to check against no variation in the buffer
                % at all.
                dvm = median(buf(ii,1));
                % mad works poorly with diff(Eh) because most
                % values are 0 so almost always the mad is 0.  This
                % is known as "implosion" of the scale estimate s and happens when at least 50% of
                % the data has the same value.  It is common for
                % quantized data.  Since Eh is quantized, diff(Eh)
                % often yields zero values.  The Eh signal itself
                % is generally free of outliers.  We're using this
                % here to decide whether the real data is abnormal
                % relative to background.  So a better estimate of
                % the first difference should be fairly easy to
                % obtain.
                s = 1.4826*median(abs(buf(ii,1)-dvm));
                a = abs(d(1)-dvm)./s > ...
                    norminv(1-PF/2,0,1);

                if (a)
                    disp('outlier');
                    b = 2^4;
                end
                
            end                
            function b = mask_depth_max(obj,d,buf,MAX);
                b = 0;
                if d(5) > MAX
                    b = 2^5;
                end
            end
            function b = mask_depth_minmax(obj,d,buf,MIN,MAX);
            % convenience function for likely common combined
            % filter to avoid the need for two stacks.
                b = obj.mask_depth_min(d,buf,MIN);
                b = b + obj.mask_depth_max(d,buf,MAX);
            end
            function b = mask_time_max(obj,d,dbug,MAX);
                b = 0;
                if d(2) > MAX
                    b = 2^6;
                end
            end


            % All filters will act only on unmasked values.
            % Should I force this or implement it within each
            % filter?   Probably implement in each filter.  Maybe
            % some filters will want to use (certain) masked
            % values.  There should be an easy function though to
            % return the unmasked values.  Filters that assume
            % regular sampling will have a hard time with some
            % masked values.
            
            function [d,t] = get_unmasked_data(obj)
            % This comes out in oldest:youngest order.

                ib = obj.buffer.fst:obj.buffer.lst;
                mask = obj.buffer.raw(ib,7);
                ii = ~logical(mask);
                d = obj.buffer.raw(ib,1);
                d = d(ii);
                t = obj.buffer.raw(ib,2);  % assumed indecies.
                t = t(ii);
            end
            function [d] = get_last_data(obj)
                d = obj.buffer.raw(obj.buffer.lst,:);
            end
            function [dbuf] = get_unmasked_dbuf(obj)
            % this is a more useful version of the above.
            % Retaining the above to avoid breaking previously
            % implemented prototype accumulators.
                
                dbuf = obj.buffer.raw(obj.buffer.fst: ...
                                      obj.buffer.lst,:);
                ii = ~logical(dbuf(:,7));
                dbuf = dbuf(ii,:);
            end
            
            %% Filters %%%%%%%%%%%%%%%%
            % 
            % filters manipulate data rather than
            % just marking it.  Some masking operations in
            % preparation for accumulation are simpler as filters -
            % e.g. NaN data that should be ignored by accumulators.
            % Or just add to the accumulators a mask to check
            % against? No because sometimes we'll want to xor, not,
            % and, etc operations, and that will be easier in the
            % filter step.
            
            function f = filter_detrend_cov(obj,d,dbuf)
            % detrend wrt covariate rather than time.  what should
            % this be called?  "trend" is reserved for time. defer.
            % this is removing a background profile.
                
            % these parameters will not be part of DslTsStack, but
            % rather part of a filter class with its set_parameters
            % method used to set these parameters.
                b = obj.trend_intercept;
                m = obj.trend_slope;

                d = obj.get_last_data();  
                if logical(d(7)) % masked
                    f = NaN;
                    return
                end
                
                f = d(1) - (d(8)*m + b);
                
            end
            function f = filter_mask_not(obj,d,dbuf,mask)
            % NaN-out all data that isn't masked by exactly mask.
                if (d(7) == mask)
                    f = d(1);
                else
                    f = NaN;
                end
            end
            function f = filter_median(obj,d,dbuf,len)
            % median filter (not rejector - that would have to be
            % implemented as a mask)
                if (len > obj.buflen)
                    % how should this be handled for real?  It
                    % can't be a compile-time fail because len will
                    % be an inifile thing, presumably?
                    warning(['buflen is insufficient to apply desired ' ...
                             'median filter']);
                    return
                end
                [d,t] = obj.get_unmasked_data();
                if (length(d) < len)
                    disp(['not enough unmasked points yet in ' ...
                          'buffer']);
                    f = NaN;
                    return
                end
                % causal filter so there is a phase shift.
                % unclear whether should fail on points that are
                % old or what (since masked values may occur anywhere).
                dd = d(end:-1:end-len+1);  % d is arranged oldest first.
                f = median(dd);
            end
            
            % Example filter - numerical derivative - simplest
            % possible form (not for actual use).
            % Using this at minimum requires adding noise to deal
            % with quantization.
            function f = filter_diff(obj,d,dbuf)

            % get unmasked data.
                [d,t] = obj.get_unmasked_data();
                
                % check for buffer being sufficiently full.
                % assuming scalar data is:
                % (data,time,northing,easting,depth,...)
                f = NaN;
                if length(d) > 1
                    dt = t(end)-t(end-1);
                    if  dt > 0
                        f = (d(end)-d(end-1))/dt;
                    end
                end
            end

            % derivative estimated from polynomial fit to latest n
            % samples in buffer.  Slow implementation of
            % Savitsky-Golay filter but tolerant of unevenly spaced
            % data.
            % M: order of polynominal fit
            % N: number of samples to apply
            %
            % This introduces phase.
            %
            % Seems to only work well for a line fit, otherwise
            % polyfit complains about bad scaling of data, which
            % could be corrected, maybe.  Using this to deal with
            % quantization noise is probably not a great idea.
            function f = filter_diff_poly(obj,d,dbuf,M,N)

            [dbuf] = obj.get_unmasked_dbuf();

                if size(dbuf,1) < N
                    disp('Not enough samples for polynomial fit');
                    f = NaN;
                    return
                end
                p = polyfit(dbuf(end:-1:end-N+1,2),dbuf(end:-1:end-N+1,1),M);

                % polynomial: 
                %   P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1)
                % derivative of polynomial:
                %   P(1)*X^(N-1)*N + P(2)*X^(N-2) + .... P(N-1)*x + P(N)
                % e.g. for M=3:
                %      x = p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4)
                %   dxdt = p(1)*x^2*3 + p(2)*x*2 + p(3)
                % for M=1 this should just give the slope of the line-fit.
                x = dbuf(end-floor(N/2),2);  % center value (phase shift).
                dxdt = 0;
                %p
                %x
                for mm = 1:M
                    dxdt = dxdt + p(mm)*x^(M-mm)*(M-mm+1);
                    % M = 3:
                    % mm =1       p(1)*x^(2)*3 
                    % mm =2       p(2)*x^(1)*2
                    % mm =3       p(3)*x^(0)*1 = p(3)
                end
                f = dxdt;
            end
            
            % convenient helper function for accumulators.
            % Ignores NaNs to facilitate use of filtering step as
            % final mask.  
            function [d,t] = get_filtered_data(obj)

                dd = obj.buffer.raw(obj.buffer.fst: ...
                                     obj.buffer.lst,:);
                ii = ~isnan(dd(:,1));
                dd = dd(ii,:);
                if size(dd,1) < 1
                    d = [];
                    t = [];
                    return
                end
                d = dd(:,1);
                t = dd(:,2);
                % really we're often going to want all the
                % ancillary data as well.  Just prototyping here.
            end
            function [dbuf] = get_filtered_dbuf(obj)
            % this is a more useful version of the above.
            % Retaining the above to avoid breaking previously
            % implemented prototype accumulators.
                
                dbuf = obj.buffer.raw(obj.buffer.fst: ...
                                      obj.buffer.lst,:);
                ii = ~isnan(dbuf(:,1));
                dbuf = dbuf(ii,:);
            end
            function [d] = get_filtered_last(obj)
                d = obj.buffer.raw(obj.buffer.lst,:);
                if isnan(d(9)) % ignore if nan'd.
                    d = [];
                end
            end

            
            %% Accumulators %%%%%%%%%%%%%%%%%%
            %
            % To facilitate use of filters as final masks,
            % accumulators should ignore NaNs.
            % nsamp should be an input parameter but just use
            % acclen for this prototype.
            function acc = get_accumulated_dbuf(obj)
                ii = ~isnan(obj.accumulated(:,1));
                acc = obj.accumulated(ii,:);
            end
            function accumulate_sample_algR(obj)
                nsamp = obj.acclen;
                [d] = obj.get_filtered_last();
                if isempty(d)
                    disp('last filtered sample invalid');
                    return;
                end
                if obj.num_accumulated < nsamp
                    % then definitely add this sample.
                    disp(['nsamp > number of filtered samples available. ' ...
                          ' Adding this data.']);
                    obj.num_accumulated = obj.num_accumulated+1;
                    obj.accumulated(obj.num_accumulated,:) = d;
                    obj.num_samples = obj.num_samples + 1;
                    return
                end
                % otherwise should we replace a sample?  Probably
                % this should be weighted by time so that samples
                % end up being distributed somewhat evenly over
                % time.  Eventually this will still replace all old
                % samples so this is probably a pretty bad
                % iterative sampler.   But turns out this is a real
                % method and has some mathematical guarantees.  It
                % is the most common "Resevoir sampler" and is
                % called Algorithm R by Jeffrey Vitter.
                obj.num_samples = obj.num_samples + 1;
                r = randi(obj.num_samples,1);  % r \in [1:obj.num_samples]
                if r <= obj.acclen
                    % then replace
                    obj.accumulated(r,:) = d;
                end
                
                
            end
            
            function accumulate_lsline(obj)
            % least squares line fit to data vs. covariate.
            % trend is vs. time so this should be called something
            % else...  _lsline.  
            % a robust alternative would be to use Theil Sen
            % algorithm.
                
            % generate (and store/update) a random sample.
                obj.accumulate_sample_algR();
                % The above is the only part that is actually an
                % accumulator.  Below is just a utility function
                % that should be called outside this class.
                %
                % It should still be implemented as part of the
                % library though because users shouldn't have to
                % deal with things like an incompletely filled
                % accumulator.  They don't - just use the
                % get_accumulated method.
                %
                % leaving this here for now because of time.
                
                % LS fit a line.
                d = obj.get_accumulated_dbuf();
                if size(d,1) < 2
                    disp('not enough samples to fit a line');
                    return
                end
                % covariates could be nan if they went in bad in
                % which case this won't work so probably should
                % ignore all data that enters with any element
                % being nan.
                A = [d(:,8) ones(size(d,1),1)];
                b = d(:,1);
                x = A\b;
                obj.trend_intercept = x(2);
                obj.trend_slope = x(1);
                
                % also compute prediction interval.
		%s = std(
                
                
            end
            
            function accumulate_tail_dt(obj,DT)
            % best obj.acclen values, with a minimum of DT temporal
            % separation.

                a = obj.accumulated;
                
                d = obj.get_filtered_last();
                if isempty(d)
                    disp('empty');
                   return
                end
                
                iw = [];
                for n = 1:obj.acclen
                    if isnan(a(n,1));
                        % then just stick this value in here and
                        % return.
                        obj.accumulated(n,:) = d;
                        disp(sprintf('filling %d/%d...',n,obj.acclen));
                        return
                    end
                    if d(2) >= abs(a(n,2))-DT ...
                            & d(2) <= abs(a(n,2)+DT)
                        % then it's in the window of an already
                        % saved outlier.  Check if this one is more
                        % intense and replace it if so.
                        %
                        % That's not enough.  It may have multiple
                        % outliers in its window.  Desired behavior
                        % isn't clear in that case.  Replace the
                        % weakest one?
                        % 
                        % A smarter plan may be to have a rather
                        % larger buffer here and follow this with a
                        % chain that does only temporal clustering.
                        % Filtering differently may also fix this.
                        disp('found other in window');

                        iw = [iw n];
                    end
                end
                
                if ~isempty(iw)
                    % if bigger than the smallest in the window, replace it
                    [dm,im] = min(abs(a(iw,1)));
                    ii = iw(im);
                    if abs(d(1)) > dm
                        obj.accumulated(ii,:) = d
                        disp('replaced in window');
                    end
                else 
                    % If a previously accumulated anomaly has been
                    % replaced then there is nothing more to do.  If
                    % not, then then any new ones have to get compared
                    % against all anomalies thus far acquired.
                    [dm,im] = min(abs(a(:,1)));
                    if abs(d(1)) > dm
                        obj.accumulated(im,:) = d;
                        disp('replaced global min.');
                    end
                end
                
            end                
            
        end

    end
    
