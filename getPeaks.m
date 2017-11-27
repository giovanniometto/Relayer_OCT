function [V, L] =  getPeaks(P, Options)

% GETPEAKS returns the value and location of the highest 'nFirst' peaks in 
% profile P sorted by value or by location if the option sortByLocation is
% 'true'. Peaks are selected from those higher than a treshold value
% defined by the ratio with the highes peak ([0 ... 1]).
%
% inputs,
%   P : the profile array of values
%   Options : A struct with all snake options
%
% outputs,
%   V : values of the peaks
%   L : location of the peaks as the index in the array
%
% options (general),
%  Option.nFirst : number of the first peaks wanted, '0' for all peaks (default)
%  Options.sortByLocation : if true, sort by location, default by value
%  Options.Order : 1 ascending / -1 descending (default)
%  Options.Treshold : only peaks higher than Treshold * value of the
%                     highest peak in the profile, default no treshold
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision:

V = [];
L = [];

defaultoptions=struct('nFirst',0,'returnSortedByLocation',false,'Order',-1,'Treshold',0);
if(~exist('Options','var'))
    Options = defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options)))
        warning('snake:unknownoption','unknown options found');
    end
end

% find max peaks 
[ maximaPks , maximaLocs ] = findpeaks(P);
Pks = [ maximaPks maximaLocs ];
if isempty(Pks), return, end

% treshold peaks
if Options.Treshold ~= 0
    minValue = Options.Treshold * max(maximaPks(:));
    Pks( Pks(:,1) < minValue, : ) = [];
end

% sort by highest value
sortedPks = sortrows(Pks , -1 );

% get highest nFirst
if Options.nFirst ~= 0 && Options.nFirst < size(sortedPks,1)
    sortedPks = sortedPks(1:Options.nFirst,:);
end

% return sorted (by value or location)
sortedPks = sortrows(sortedPks , Options.Order * (1+(1*Options.returnSortedByLocation)) );
V = sortedPks(:,1);
L = sortedPks(:,2);
