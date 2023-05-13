% Emulate processing of hydrothermal data documented in jakuba07
% for dive abe126.  Processing is not exactly identical because
% some of the original processing was batch.
%
% Revision History
% 2016-11-28    mvj    Created.

MOVIE = false;

% load data.  Struct htp contains raw and processed sensor data.
load('/data/abe2004/abe126/abe126.mat','htp','nav','state','sci');

% Eh: 
% 1) segment out near-surface data
% 2) numerical derivative
% 3) windowed-hampel identifier for negative outliers.
% 4) accumulate tail based on magnitude of negative slope.


% Implementation in terms of stacks:
% A. (1), (2)
% B. (3) - acts on masked dEhdt; (4) does temporally separated tail
% of anomalies.

e1 = dehdt_stack();
e2 = eh_stack();
e3 = ehanom_stack();

% htp is clipped in time for on-bottom, so use the more raw version
% in sci.

T = 56700;  % ignore ascent.
            % T = 20000;  % just first part.
t0 = sci.t-sci.t(1);

t1 = NaN*ones(T,1);
v1 = NaN*ones(T,1);
f1 = NaN*ones(T,1);
m1 = NaN*ones(T,1);
t2 = NaN*ones(T,1);
v2 = NaN*ones(T,1);
f2 = NaN*ones(T,1);
m2 = NaN*ones(T,1);
t3 = NaN*ones(T,1);
v3 = NaN*ones(T,1);
f3 = NaN*ones(T,1);
m3 = NaN*ones(T,1);

% interpolate onto data timebase.
northing = interp1(nav.t,nav.y,sci.t);
easting = interp1(nav.t,nav.x,sci.t);
depth = interp1(nav.t,nav.z,sci.t);
altitude = interp1(state.t,state.altitude,sci.t);

% input mask
mask = false(size(northing));

% input covariate(s?).  
covariate = NaN*ones(size(northing));

if MOVIE
    video = VideoWriter('eh.avi');
    open(video);
end


for n=1:T

    % progress
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d/%d',n,T);
    
    % run processing chain.
    d1 = [sci.eh(n), t0(n), northing(n), easting(n), depth(n), ...
          altitude(n), mask(n), covariate(n)];
    e1.add_data(d1);
    d2 = e1.get_last_data();  % this will already have a bunch of
                              % masks applied.  Ok, I guess.
    d2(1) = d2(9); % filtered value is input into next chain.
    d2(9) = [];  % awkward.
                 % filtered value may not yet be valid, e.g. diff
                 % needs two values to be in the buffer, or it may
                 % have failed for other reasons.  Probably the
                 % class above should mark this.  Do it here for
                 % now.
    if isnan(d2(1))
       d2(7) = 2^0;  % use input mask to indicate value is bad. 
    end
    e2.add_data(d2);
    d3 = e2.get_last_data();
    d3(1) = d3(9); % filtered value is input into next chain.
    d3(9) = [];
    e3.add_data(d3);
    
    % save results for plotting purposes.
    t1(n) = e1.buffer.raw(e1.buffer.lst,2);
    v1(n) = e1.buffer.raw(e1.buffer.lst,1);
    f1(n) = e1.buffer.raw(e1.buffer.lst,9);
    m1(n) = e1.buffer.raw(e1.buffer.lst,7);
    t2(n) = e2.buffer.raw(e2.buffer.lst,2);
    v2(n) = e2.buffer.raw(e2.buffer.lst,1);
    f2(n) = e2.buffer.raw(e2.buffer.lst,9);
    m2(n) = e2.buffer.raw(e2.buffer.lst,7);
    t3(n) = e3.buffer.raw(e3.buffer.lst,2);
    v3(n) = e3.buffer.raw(e3.buffer.lst,1);
    f3(n) = e3.buffer.raw(e3.buffer.lst,9);
    m3(n) = e3.buffer.raw(e3.buffer.lst,7);

    
    if ~mod(n,100) && n > 5000
        it = 1:n;
        figure(1); clf reset
        subplot(211)
        plot(sci.t(1)+t1(it),v1(it)*1000);
        hold on;
        ii = m2==2^4;
        oeh = interp1(t1(it),v1(it),t3(ii));
        plot(sci.t(1)+t3(ii(it)),oeh*1000,'.');
        a = e3.get_accumulated_dbuf;
        % anomaly is in dEhdt not Eh, so find the hits.
        aeh = interp1(t1(it),v1(it),a(:,2));
        plot(sci.t(1)+a(:,2),aeh*1000,'ko');
        grid on;
        legend('Raw','Anomalous','Most Intense');
        ylabel('Eh (ORP) [mV]');
        
        subplot(212)
        plot(sci.t(1)+t1(it),f1(it)*1000,'-')  
        hold on
        ii = m3==2^4;
        plot(sci.t(1)+t3(ii(it)),v3(ii(it))*1000,'.')
        plot(sci.t(1)+e3.accumulated(:,2),e3.accumulated(:,1)*1000,'ko')
        grid on;
        %legend('Raw','Anomalous','Most Intense');
        ylabel('dEh/dt [mV/s]');
        twin([1.0e+09 *1.095272616535433    sci.t(end)]);
        tticklabel('abs',2*3600);
        
        drawnow
        
        if MOVIE 
            writeVideo(video,getframe(gcf));
        end
        
    end
    
    
end

if MOVIE
    close(video);
end

return

figure(1); clf reset
subplot(211)
plot(sci.t(1)+t1,v1*1000);
hold on;
ii = m2==2^4;
oeh = interp1(t1,v1,t3(ii));
plot(sci.t(1)+t3(ii),oeh*1000,'.');
a = e3.get_accumulated_dbuf;
% anomaly is in dEhdt not Eh, so find the hits.
aeh = interp1(t1,v1,a(:,2));
plot(sci.t(1)+a(:,2),aeh*1000,'ko');
grid on;
legend('Raw','Anomalous','Most Intense');
ylabel('Eh (ORP) [mV]');

subplot(212)
plot(sci.t(1)+t1,f1*1000,'-')  
hold on
ii = m3==2^4;
plot(sci.t(1)+t3(ii),v3(ii)*1000,'.')
plot(sci.t(1)+e3.accumulated(:,2),e3.accumulated(:,1)*1000,'ko')
grid on;
%legend('Raw','Anomalous','Most Intense');
tticklabel('abs',3600);
ylabel('dEh/dt [mV/s]');
twin(1.0e+09 *[   1.095272616535433    1.095289907874016]);
