function [xarr,plotarr] = postproc(obj,varargin)
% plots full branches
% varargin{1} should be bcarr
% varargin{2} can specify a preference between 'E' 'D' and 'I'

close all

switch length(varargin)
    case 0
        bcarr = obj;
        plotvar = 'E';
    case 1
        bcarr = varargin{1};
        plotvar = 'E';
    case 2
        bcarr = varargin{1};
        plotvar = varargin{2};
end

% direction.names in BC determine the possible plotting directions
% (more or less) -- so make sure you are calling this from the
% right BC
Warr = [bcarr.Wsave]; 
xarr = zeros(length(obj.direction.names),length(Warr));
for j = 1:length(obj.direction.names)
    xarr(j,:) = [Warr.(obj.direction.names{j})];
end
inds = (max(xarr,[],2)-min(xarr,[],2)) ~= 0;

if isempty(inds)
    disp('There is nothing to plot here, boss.')
    return
else
xarr = xarr(inds,:);
plotbc = [bcarr.to_plot];
plotarr = [plotbc.(plotvar)];

names = plotnames(obj.direction.names);
names = names(inds);
if length(names) < 2
    figure
    plot(xarr,plotarr,'.-k','linewidth',1,'markersize',11);
    xlabel(names{1})
    ylabel(plotvar)
else
    figure
    plot3(xarr(1,:),xarr(2,:),plotarr,...
        '.-k','linewidth',1,'markersize',11);
    xlabel(names{1})
    ylabel(names{2})
    zlabel(plotvar)
    grid on
end
end

end

function names = plotnames(dirnames)
parnames = {'l' 'beta' 'Re' 'Wi' 'Sc' 'epsilon' 'lambda' 'k'};
displaynames = {'$L_{max}$' '$\beta$' '$Re$' '$Wi$' '$Sc$'...
    '$\varepsilon$' '$\lambda$' '$k$'};
names = cell(size(dirnames));
for j = 1:length(dirnames)
ind = find(strcmp(parnames,dirnames{j}));
names{j} = displaynames{ind};
end
end
