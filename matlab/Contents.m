% Matlab prototype of basic timeseries processing for 
% acoustic telemetery of chemical data, especially anomalies
%
% Two processing chains were implemented.  These approximate
% the batch processing applied to certain historical hydrothermal 
% data collected by ABE:
%
%   * Eh/ORP anomalies indicative of buoyant hydrothermal water
%   * OBS significantly above background to indicate non-buoyant 
%     hydrothermal plume water.
%
% These two examples are representative of two broad class of 
% anomalies: those associated with point features, and those 
% that occur over a relatively large portion of the survey 
% domain.
%
% Quick start:
% ------------
% Data necessary to execute the chains are available here:
% @@@
% Download the this data to the directory containing this file.
%
% To execute the Eh/ORP processing chain at the matlab prompt 
% >> abe126_eh
%
% To execute the OBS processing chain
% >> abe126_obs
%
% Data format:
% ------------
% Input data consists of scalar sensor data plus covariates:
%
% [raw,time,northing,easting,depth,altitude,mask,covariate]
%
% Time and position facilitate operations such as temporal 
% and spacial clustering.  An input mask can be used to, 
% for example, ignore data from a descent.  The generic 
% covariate may be, for example, potential density to allow 
% removal of background variation.
% 
% Comments on implementation:
% ---------------------------
% The implementation here is a prototype.  Although the 
% intention was to write the code in such a way as to 
% facilitate implementation in C++, there are some stuctural
% problems that compromise usability.  At a minimum
% the collection of low-level processing functions should
% be useful.
%
% The processing chains are implemented as classes derived from 
% DslTsStack.  DslTsStack requires users to implement three
% functions, apply_mask, apply_filter, and accumulate.  These 
% three functions appear to be adequate for basic timeseries 
% processing:
%   apply_mask - tag data, for example as anomalous or invalid
%   apply_filter - modify the raw data, e.g. differentiate
%   accumulate - compute long-term statistics (means, trends),
%                retain samples, or other operations that 
%                depend on more than an a short history of
%                recent data.
%
% DslTsStack contains a circular buffer of a size specified
% at instantiation to facilitate operations on recent data.
% 
% DslTsStack implements all processing functions (masks, 
% filters, accumulators).  Derived classes implement the
% operations by specifing a particular method from DslTsStack
% to apply for each of the steps above.  Wrapper scripts
% chain together multiple such derived stacks into complete
% processing chains.
%
% There are at least two major issues with the prototype
% implementation here:
%  (1) Derived classes are required for each stack, limiting
%      reusability
%  (2) Each stack must consist of exactly one mask, filter
%      and accumulator
%
% An important architectural improvement would be to 
% allow direct use of masks, filters, and accumulators.  
% This would make the entry of variable length lists of 
% parameters or evolving parameters happen in user code, 
% i.e. at the scripts level here, eliminate the awkward
% use of class data to pass parameters between stacks, and
% eliminate the need for user-implemented derived classes.  
% The circular buffer should probably remain encapsulated 
% in a class that abstracts insertion and retrival of data.  
%
% Processing functions (implemented in DslTsStack.m):
% --------------------
% 
%   Masks:
%     min/max value, depth, and time
%     Hampel identifier (outlier detection)
%   Filters:
%     remove background (trend wrt a covariate)
%     median filter
%     first order numerical derivative
%     differentiate from a polynomial fit
%   Accumulators
%     resevoir sample (algorithm R)
%     least squares line fit (to resevoir sample)
%     temporally separated extrema 
%   
% 
% Matlab functions:
% -----------------
%   circVbuf    - circular buffer (from Matlab Exchange)
%
%   dslTsStack  - top level class encapsulating data
%                 and processing functions.
%
%   abe126_eh            - Eh processing script, ABE dive 126.
%   dehdt_stack          - Differentiates Eh
%   eh_stack.m           - Identifies short-term outliers
%   ehanom_stack.m       - Retains most significant temporally-
%                          separated outliers
%   abe126_obs           - OBS processing script, ABE dive 126.
%   obs_stack.m          - de-spike (median filter)
%   obsprofile_stack.m   - compute profile wrt potential density
%   obsanom_stack.m      - compute anomaly relative to background
%


