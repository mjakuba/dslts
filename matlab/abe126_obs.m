% Emulate processing of hydrothermal data documented in jakuba07
% for dive abe126.  Processing is not exactly identical because
% some of the original processing was batch.
%
% This is for OBS, which is slightly interesting because NBP anomaly is
% defined relative to a trend and using prediction interval.
%
% Revision History
% 2016-11-28    mvj    Created.


% load data.  Struct htp contains raw and processed sensor data.
load('/data/abe2004/abe126/abe126.mat','star','nav','sci','htp','state');

% Background pdens:
% 0) (done externally) compute pdens from in situ C,T,P,latitude.
% 1) segment by depth (to some nominal height above expected plume
%    height - looks like 1625 was used for abe126).  This is needed
%    to use sigma as a covariate in rest of OBS stack.  Seems
%    awkward to do it here, since normally you would need the
%    variate itself to do this segmentation.
% 2) probably no filtering required
%
%
% pdens can be processed in a single stack.
%
% Seems there is no need at all for a pdens stack.  The use of
% pdens as a covariate for OBS can use the same segmentation
% applied to the variate, either just by depth or something that
% actually attempts to identify the NBP from the data.

% OBS 
% 0) interpolate pdens onto obs timebase and add as covariate.
% 1) segment by depth
% 2) median filter to remove spikes (10-pt)
% 3) detrend (from sigma trend provided by separate stack)  Input
% of changing variables.  Not sure how to handle this...
% 4) compute prediction interval 
% 5) apply it prediction interval to NBP OBS classification.


% Implementation in terms of stacks:
% A. (1), (2)
% B. (3), (4) filter only (detrend)
% C. (5)
% of anomalies.

e1 = obs_stack();  % this should be obspre_stack
e2 = obsprofile_stack();
e3 = obsanom_stack();

T1 = 19000; % descent plus some of the survey including the BP hit.
t0 = sci.t-sci.t(1);

t1 = NaN*ones(T1,1);
v1 = NaN*ones(T1,1);
f1 = NaN*ones(T1,1);
m1 = NaN*ones(T1,1);
c1 = NaN*ones(T1,1);
t2 = NaN*ones(T1,1);
v2 = NaN*ones(T1,1);
f2 = NaN*ones(T1,1);
m2 = NaN*ones(T1,1);
c2 = NaN*ones(T1,1);
t3 = NaN*ones(T1,1);
v3 = NaN*ones(T1,1);
c3 = NaN*ones(T1,1);
f3 = NaN*ones(T1,1);
m3 = NaN*ones(T1,1);


for n=1:T1

    % temporarily subsample for speed.
    %if mod(n,10)
    %    continue;
    %end
    % progress
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d/%d',n,T1);
    
    
    % interpolate onto data timebase.
    northing = interp1(nav.t-nav.t(1),nav.y,t0(n));
    easting = interp1(nav.t-nav.t(1),nav.x,t0(n));
    depth = interp1(nav.t-nav.t(1),nav.z,t0(n));
    altitude = interp1(state.t-state.t(1),state.altitude,t0(n));
    
    % interpolate sigma1 onto data timebase.  Probably we need a
    % utility function for this, or in practise will probably just
    % take the last value - that will happen external to dslts though.
    sigma = interp1(star.t-star.t(1),star.sigma1,t0(n));
    
    % input mask
    mask = false;
    
    % input covariate(s?).  
    covariate = sigma;
    
    % run processing chain.
    
    % 1. mask surface data and de-spike.
    d1 = [sci.obs(n), t0(n), northing, easting, depth, ...
          altitude, mask, covariate];
    if (sum(isnan(d1)) > 0) 
        continue;
    end    
    e1.add_data(d1);
    
    % 2. estimate background profile and prediction interval.
    % this assumes there is only one covariate.
    dd = e1.get_filtered_last();
    if isempty(dd)
        % filtered output is masked.
        continue;
    end
    d2 = dd(1:8); % maybe just retain the filtered value to avoid
                  % this extra line?  Or just make a helper
                  % function to do this for you in one step since
                  % sending the filtered value into the next step
                  % is going to be pretty common.
    d2(1) = dd(9); % use filtered value as input to next stage.
    e2.add_data(d2);
    % library will have utility functions for computing and LS line
    % and prediction interval, and detrending data.  But these will
    % be called externally in implementations.
    % already implemented LS line-fit as an accumulator so just
    % going to use that here.  But implement detrend and pred
    % interval here so that it's closer to how this will actually
    % work.
    aa2 = e2.get_accumulated_dbuf();  % nan-free.
    m = e2.trend_slope; % computed outside in future
    b = e2.trend_intercept; % computed outside in future 
    if isnan(m)
        continue
    end
    sigma = aa2(:,8);
    obsf = aa2(:,1);
    obsfd = obsf - (m*sigma + b); % remove background profile.
    % @@@ pretty sure the above is making things worse.
    sv = std(obsfd);  % this is all that is needed for pred
                      % interval vs sigma.
    
    % now there is a problem.  How does the trend-line computed
    % in the step above get into the detrender?  Easy enough to get
    % out of e2.  But these are effectively parameters that vary
    % over time.  The accumulator should be general.  But there
    % needs to be an interface that says "use these parameters for
    % the accumulator when adding data."  So it needs to be part of
    % add_data.  Or perhaps a separate function that has to be
    % called beforehand to update the parameters.  remember: things
    % like obs_stack we want to not require (i.e. the user should
    % only rarely need to derive classes from DslTsStack.  But the whole stack
    % will always be compiled by the user, so this additional step
    % doesn't break the intended usage so long as each DslTsStack
    % is limited to one mask, filter, and accumulator.
    % 
    % The way to get parameters into everything is to have each
    % mask, filter, and accumulator define a set_parameters method.

    % 3. Detrend to compute the OBS anomaly.  This is a
    % just a filter parameterized by slope and intercept estimated
    % from the output of (2).  Except this is applied to all data,
    % except descent.
    % setting the trend parameters would normally be done by
    % calling the set parameters method of the accumulator in e3,
    % but accumulator is not implemented as a class in this prototype.
    e3.trend_slope = m;
    e3.trend_intercept = b;
    e3.add_data(d2);  % yes, d2.  The stack of (2) does not produce
                      % data we use subsequently.
    % note - the filtered output is detrended wrt sigma, not time.
    % Right... so what is the anomaly vs. time in native OBS units?  That's what we
    % want in the filtered output.
    
    
    
    % 4. find NBP from OBS data.  (Compute the prediction
    % interval).  This again operates on a sample of data from
    % a (known) depth interval above the NBP.  Maybe it should be
    % combined with step 2, or use the accumulated data from step 2
    % to insure the statistics are consistent.
    %
    % Defer - the output is not clear and low priority.  NBP tends
    % to cover most of survey and so doesn't lend itself to
    % incremental statistics.  The way this would be handled is to
    % send up the prediction interval and the detrended OBS data
    % from which a topside plotter could classify the data.  
    %
    % Right - so the way to do this is to apply the prediction
    % interval as a mask here and send that up.  Topside can still
    % potentially do better because it has as much data as it has
    % collected, whereas this can't reprocess old data.
    
    
    % 5. find BP from OBS data.
    
    
    % save results for plotting purposes.
    t1(n) = e1.buffer.raw(e1.buffer.lst,2);
    v1(n) = e1.buffer.raw(e1.buffer.lst,1);
    f1(n) = e1.buffer.raw(e1.buffer.lst,9);
    m1(n) = e1.buffer.raw(e1.buffer.lst,7);
    c1(n) = e1.buffer.raw(e1.buffer.lst,8);
    t2(n) = e2.buffer.raw(e2.buffer.lst,2);
    v2(n) = e2.buffer.raw(e2.buffer.lst,1);
    f2(n) = e2.buffer.raw(e2.buffer.lst,9);
    c2(n) = e2.buffer.raw(e2.buffer.lst,8);
    m2(n) = e2.buffer.raw(e2.buffer.lst,7);
    t3(n) = e3.buffer.raw(e3.buffer.lst,2);
    v3(n) = e3.buffer.raw(e3.buffer.lst,1);
    f3(n) = e3.buffer.raw(e3.buffer.lst,9);
    m3(n) = e3.buffer.raw(e3.buffer.lst,7);
    c3(n) = e3.buffer.raw(e3.buffer.lst,8);
    
end

figure(1); clf reset;
plot(t1,v1,'.')
hold on
plot(t2,f2,'o')
legend('Raw OBS','Descent');
xlabel('Time (s)');
ylabel('OBS mV');

figure(2); clf reset;
plot(c1,v1,'o')
hold on
plot(c2,v2,'.')
plot(c3,f3,'x')
legend('Raw','Med. Filtered','Background Rm''d');
