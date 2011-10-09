
function RGBmap = colorGray(numLevels,debugplot)

% BRIEF: creates a rainbor color colormap which also looks good in greyscale
%
% SYNTAX:
%     RGBmap = colorGray(numLevels)
%
%OUTPUTS
%
%INTPUTS
%    numlevels - number of levels in the output colormap
%    plot - boolean, if true outputs a plot showing the range of colors
%    produced and lineplot of the range of grayscale produced verifying
%    grayscale linearity.
%
%****************************************************
%Author: Alexandre R. Tumlinson, October 23, 2006
%****************************************************
%REVISIONS
%16, Nov. 2006 Peder Axensten - sent great suggestions and some nicely rewriten code.
% make the numLevels input optional. Also optimized by vectorizing
% code.  Added the debug plotting options.
%****************************************************

if(nargin < 1), numLevels= size(colormap, 1); end	% Same number of  colors as the present color map. 
if(nargin < 2), debugplot=0; end	% turn off debug plotting by default

%first make a Jet map and then scale it so that it will look good in gray
lightest=0.05;
% Some constants to trim the greymap so lines arent too light.
% Also limit dark end because color sarted looking weird.
offsetbrightlimit=	floor(numLevels*lightest);
offsetdarklimit=	ceil(numLevels*0.25);
grayInt=			gray(numLevels+offsetbrightlimit+offsetdarklimit);
grayInt=			grayInt(offsetdarklimit+1:end-offsetbrightlimit,:);
grayInt=			sum(grayInt,2)/3.4;

% Do something similar for the color map
offsetbluelimit=	floor(numLevels*0.1);
offsetredlimit=		ceil(numLevels*0.0);
jetmap=				jet(numLevels+offsetbluelimit+offsetredlimit);
jetmap=				jetmap(offsetbluelimit+1:end-offsetredlimit,:);
jetmapInt=			jetmap*[0.2989; 0.5870; 0.1140]; %the gray level of the jetmap

scalefactor=		grayInt./jetmapInt;
RGBmap=				jetmap.*repmat(scalefactor, 1, 3);

%find elements greater than 1 and distribute to other colors
%attempts to maintain grey level on redistribution
%redistribution percentages are a litte arbitrary, but it seems to work
% Find elements greater than 1 and distribute to other colors,
% attempts to maintain grey level on redistribution.
% Redistribution percentages are a litte arbitrary, but it seems to work
for trials = 1:5
    isdone=				true;
    bigs=				(RGBmap(:,1) > 1);
    if (any(bigs))
        isdone=				false;
        surplus=			RGBmap(bigs,1) - 1;
        RGBmap(bigs,1)=		1;
        RGBmap(bigs,2)=		RGBmap(bigs,2) + surplus*0.85*(0.2989/0.5870);	% Goes to salmon
        RGBmap(bigs,3)=		RGBmap(bigs,3) + surplus*0.15*(0.2989/0.1140);	% Goes to salmon
        %			RGBmap(bigs,2)=		RGBmap(bigs,2) + surplus*0.00*(0.2989/0.5870);	% Goes to pink
        %			RGBmap(bigs,3)=		RGBmap(bigs,3) + surplus*1.00*(0.2989/0.1140);	% Goes to pink
    end
    bigs=				(RGBmap(:,2) > 1);
    if (any(bigs))
        isdone=				false;
        surplus=			RGBmap(bigs,2) - 1;
        RGBmap(bigs,1)=		RGBmap(bigs,1) + surplus*0.50*(0.5870/0.2989);
        RGBmap(bigs,2)=		1;
        RGBmap(bigs,3)=		RGBmap(bigs,3) + surplus*0.50*(0.5870/0.1140);
    end
    bigs=				(RGBmap(:,3) > 1);
    if (any(bigs))
        isdone=				false;
        surplus=			RGBmap(bigs,3) - 1;
        RGBmap(bigs,1)=		RGBmap(bigs,1) + surplus*0.00*(0.1140/0.2989);
        RGBmap(bigs,2)=		RGBmap(bigs,2) + surplus*1.00*(0.1140/0.5870);
        RGBmap(bigs,3)=		1;
    end
    if(isdone),			break;		end
end

%debugging plot option
if(debugplot)
    testarray=ones(2,numLevels);
    for count=1:numLevels
        testarray(:,count)=testarray(:,count)*count;
    end
    figure
    %create a line plot demonstrating the color range of the colormap
    subplot(1,2,1)
    set(gca,'ColorOrder',RGBmap)
    set(gca,'NextPlot','replaceChildren') %because default operation plot automatically sets default color order
    plot(testarray,'lineWidth',5 )
    title('Colors in colormap');
    ylabel('colormap index');
    set(gca,'Ylim',[0.5 numLevels+0.5]);
    set(gca,'Xtick',[]);

    % Make a plot of grayscale linearity for created colormap
    subplot(1,2,2)
    greyval=			RGBmap*[0.2989; 0.5870; 0.1140]; %the gray level of the output map
    plot(greyval);
    title('Range of grey used');
    xlabel('colormap index');
    ylabel('equivalent greylevel') %would be nice to put this on the right side of plot, any suggestions?
    set(gca,'Ylim',[0 1]);
    set(gca,'Xlim',[0.5 numLevels+0.5]);
end

